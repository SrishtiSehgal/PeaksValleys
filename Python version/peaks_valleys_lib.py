
import numpy as np
import math

###########
# VectorAllocateF - vecmat.c
# add_to_vector - arrpunt.c
# rectangular_data_matrix2file_with_a_header_using_row_and_column_indices
# split_name_extension - util.c
# multi_append2 - string-util.c
# line_segment - data_files.c
############

def sigma(amplitude, Gaussiancenter, x):
	#estimates the sigma value for 
	# y = Aexp(-(x-x0)^2/sigma)
	#if not possible value of -999 is returned

	try:
		return math.sqrt(((x-Gaussiancenter)**2)/math.log(amplitude+1)) #########why +1
	except Exception:
		return -999

def sigma2(f_amplitude, f_Gaussiancenter, f_x, f_y):
	#estimates the sigma value for 
	# y = Aexp(-(x-x0)^2/sigma)
	#if not possible value of -999 is returned

	f_sig = -999.
	f_denom = math.log(f_amplitude/f_y)

	if f_denom > 0:
		f_sig = abs(f_x-f_Gaussiancenter)/math.sqrt(f_denom)
	return f_sig

def deriv_i(vec_f_x, i_nx, i_i, i_error):
	'''
	Derived from the left
		EYE: the sampling interval for the abscissa is assumed = 1
	Input
		x: series (OFFSET = 1).
		nx: # of points of 'x'.
		i: x-coordinate of the point under analysis
	Output
		say: derived from the left or '0' if there was mess
		error = -1 if there was error
			 0 if the derivative was calculated well

	'''
	i_il = i_i-1
	
	if (i_i > i_nx) or (i_i < 0) or (i_il > i_nx) or (i_il < 0):
		i_error=-1
		return None, i_error
	else:
		i_error=0

	return vec_f_x[i_i] - vec_f_x[i_il], i_error

def deriv_d(vec_f_x, i_nx, i_i, i_error):
	'''
	Derived from the right
		EYE: the sampling interval for the abscissa is assumed = 1
	Input
		x: series (OFFSET = 1).
		nx: # of points of 'x'.
		i: x-coordinate of the point under analysis
	Output
		say: derived from the right or NULL if there was a mess
		error = -1 if there was error
		 0 if the derivative was calculated well

	'''
	i_il = i_i+1
	
	if (i_i > i_nx) or (i_i < 0) or (i_il > i_nx) or (i_il < 0):
		return None, -1
	else:
		i_error=0
		return vec_f_x[i_il] - vec_f_x[i_i], i_error

def deriv(vec_f_x, i_i, i_error):
	'''
	Derivative centered on point 'i'
	EYE: the sampling interval for the abscissa is assumed = 1
	delta x is assumed to be 1
	Input
	x: series
	i: x-coordinate of the point under analysis
	Output
	say: derivative with center in 'i'
	error = -1 if there was error
		 0 if the derivative was calculated well

	'''
	i_i2 = i_i+1
	i_i0 = i_i-1
	i_nx = int(vec_f_x[0])
	if (i_i2 > i_nx) or (i_i2 < 0) or (i_i0 > i_nx) or (i_i0 < 0):
		return None, -1
	else:
		i_error=0
		return (vec_f_x[i_i2] - vec_f_x[i_i0])/2., i_error

def second_deriv(vec_f_x, i_nx, i_i, i_error):
	'''
	Second derivative centered on point 'i'
		EYE: the sampling interval for the abscissa is assumed = 1
		Input
		x: series
		i: x-coordinate of the point under analysis
		Output
		say: second derivative with center in 'i'
		error = -1 if there was error
			 0 if the derivative was calculated well
	'''
	i_i2 = i_i+1
	i_i0 = i_i-1
	if (i_i2 > i_nx) or (i_i2 < 0) or (i_i0 > i_nx) or (i_i0 < 0):
		i_error=-1
		return None, i_error
	else:
		i_error=0
		return (vec_f_x[i_i2] - 2*vec_f_x[i_i] +vec_f_x[i_i0]), i_error

def peak_amplitude (vec_f_x, i_nx, i_ic, f_threshold, i_il, i_id):
	'''
	From the position of a supplied center (ic) of a
	given series (x) produces the succession of amplitudes under the peak of
	MAXIMUM up to the bottom level.
	The algorithm works each branch (left / right) separately.
	The criterion of stopping in a branch is: (any of them)
	1) the derivative is not calculable
	2) the derivative (module) falls below a given threshold
	3) the derivative changes its sign

	Input:
	x: vector with the time series (OFFSET = 1).
	nx: length of 'x'.
	ic: point corresponding to the center of the anomaly to study
	threshold: threshold for the derivative between two points "going down"
			from the position of the anomaly that is in 'ic'.

	Output:
	amplitude: height of the peak corresponding to 'ic' with respect to its
			 intersection with the baseline (ii, id)
	ii, id: abscissa that define the beginning and end of the peak
			 (base line)
	'''
	i_error = 0
	i_stopi = 0
	i_stopd = 0
	j = -1
	iizq = i_ic
	ider = i_ic

	while( not (i_stopi and i_stopd)): ##NOT ##AND>
		j = j + 1
		
		if (i_stopi==0):    # move fwd by 1 point on the left
			iizq = i_ic - j
			d, i_error= deriv_i(vec_f_x, i_nx, iizq, i_error) # /* derivative feasible and slope elbow */
			
			if( (i_error == -1 ) or (abs(d) <= f_threshold) or (d <= 0) ):
				#if error code is -1, 
				#or the derivative is below the threshold, 
				#or the derivative has changed sign...stop
				i_stopi = 1   #/* stop and keep the inside
		
		if (i_stopd==0):    # move fwd by 1 point on the right
			ider = i_ic + j
			d, i_error = deriv_d(vec_f_x, i_nx, ider, i_error) #/* derivative feasible and slope elbow */
			
			if( ( i_error == -1 ) or (abs(d) <= f_threshold)  or (d >= 0) ):
				i_stopd = 1 #stop and keep the inside
	#define the beginning and end of the peak
	i_il = iizq
	i_id = ider

	# peak height calculation
	i_xic = vec_f_x[i_ic]
	i_xil = vec_f_x[iizq]
	i_xid = vec_f_x[ider]

	#immediate points around peak are of the same value (still might be a peak though)
	if (ider == iizq):
		if i_xil==i_xid:
			print("\nFor this peak value {}, immediate constant values were detected around it".format(i_ic))
			return -999, i_il, i_id
		else:
			print("\nERROR: Derivative threshold too high (%f). Reduce it.\n", f_threshold)
			return None, None, None

	slope = (i_xid - i_xil)/(ider - iizq) 
	#x value  at the intersection with the peak */
	xicp  =  slope  * (i_ic  -  iizq)  +  i_xil
	amplitud = abs(i_xic - xicp)
	return amplitud, i_il, i_id

def peak_amplitude2( vec_f_x, i_nx, i_ic, f_threshold):
	'''
	From the position of a supplied center (ic) of a
	given series (x) produces the succession of amplitudes under the peak of
	MAXIMUM up to the bottom level.
	The algorithm works each branch (left / right) separately.
	The criterion of stopping in a branch is: (any of them)
	   1) the derivative is not calculable
	   2) the derivative (module) falls below a given threshold Y
	  the second derivative is positive (concavity up)
	   3) the derivative reverses the sign

	Input:
	x: vector with the time series (OFFSET = 1).
	nx: length of 'x'.
	ic: point corresponding to the center of the anomaly to study
	threshold: threshold for the derivative between two points "going down"
		from the position of the anomaly that is in 'ic'.

	Output:
	amplitude: height of the peak corresponding to 'ic' with respect to its
		 intersection with the baseline (ii, id)
	ii, id: abscissa that define the beginning and end of the peak
		 (base line)


	'''
	iizq = i_ic
	ider = i_ic
	stopi = 0
	stopd = 0
	j = 0
	i_error1=0
	i_error2=0

	while(stopi or stopd):
		j = j + 1
		if (stopi==0):    #move fwd by one point on the left
			iizq = iizq - 1
			d, i_error1 = deriv_i(vec_f_x, i_nx, iizq, i_error1)
			d2, i_error2 = second_deriv(vec_f_x,i_nx, ider, i_error2)

			if( ( i_error1 == -1 ) or ( i_error2 == -1 ) or ((abs(d) < f_threshold) and (d2 > 0)) or (d < 0)):
				stopi = 1   #stop and keep the inside

		if (stopd==0):    #move fwd by one point on the right
			ider = ider + 1
			d, i_error1 = deriv_d(vec_f_x, i_nx, ider, i_error1)
			d2, i_error2 = second_deriv(vec_f_x, ider, i_error2)

			if( ( i_error1 == -1 ) or ( i_error2 == -1 ) or ((abs(d) < f_threshold) and (d2 > 0)) or (d > 0)):
				stopd = 1   #stop and keep the inside

	#define the beginning and end of the peak
	i_il = iizq
	i_id = ider

	# peak height calculation
	i_xic = vec_f_x[i_ic]
	i_xil = vec_f_x[iizq]
	i_xid = vec_f_x[ider]

	if(ider == iizq):
	   print("\nERROR: Derivative threshold too high (%f). Reduce it.\n", f_threshold)
	   return None, None, None

	slope = (i_xid - i_xil)/(ider - iizq) 
	#x value  at the intersection with the peak */
	xicp  =  slope  * float((i_ic  -  iizq))  +  i_xil
	amplitud = float(abs(i_xic - xicp))
	return amplitud, i_il, i_id

def x_in_x1x2(f_x, f_x1, f_x2):
	# checks whether x is within the closed interval [x1, x2]
	# 1 = yes
	# 0 no
	if ((f_x >= f_x1) and (f_x <= f_x2)):
		return 1
	else:
		return 0

def deltay_pico(vec_f_x, i_ic, i_a, i_b):
	'''
	Given a point (ic) of a vectorf (x) where it exists
	one end located with base (a, b), returns the amplitude,
	measured from the string defined by its base to the peak.
	(OFFSET = 1).
	nx: size of 'x'.'
	'''
	
	m = (vec_f_x[i_b] - vec_f_x[i_a])/(i_b - i_a);
	xbase = m * (float(i_ic) - i_a) + vec_f_x[i_a];
	return float(abs(vec_f_x[i_ic]-xbase))

def simetria_pico(i_ic, i_a, i_b):
	'''
	Given a point (ic) where it exists
	one end located with base (a, b), returns the degree of symmetry
	of the position of the peak with respect to the base.

	Output
	r = symmetry factor of the peak in [0,1].

	r = MIN (fa, fb) / MAX (fa, fb);
	fa = (float) (ic-a);
	fb = (float) (b-ic);

	-1 = error in the values ​​of (ic, a, b)

   EYE: the one that (a) is the left end of the base must be respected
   and (b) the right and that (ic) is among them.
   '''

	if( ((i_ic - i_a) <= 0) or ((i_b - i_ic) <= 0)):
		return -1
	else:	
		fa = float((i_ic-i_a))
		fb = float((i_b-i_ic))
		return min(fa,fb)/max(fa,fb)

def simetria_valores(vec_f_x, i_ic, i_a, i_b):
	'''
	Given a vector (x) and point (ic) where it exists
	one end located with base (a, b), returns the degree of symmetry
	of the values ​​of the function in the position of the peak with respect to
	those of the ends of the interval that defines the base.
	(Arbitrary OFFSET, BUT NOT CHECK ANY OF THE SUBINDICES)).

	Output
	r = symmetry factor of the peak values ​​in [0,1].

	r = MIN (fa, fb) / MAX (fa, fb);
	fa = (float) fabs (x [ic] - x [a]);
	fb = (float) fabs (x [ic] - x [b]);

	-1 = error in the values ​​of (ic, a, b)

   EYE: the one that (a) is the left end of the base must be respected
   and (b) the right and that (ic) is among them
	'''

	if( (i_ic <= i_a) or (i_b <= i_ic) ):
		return -1
	else:
		fa = float(abs(vec_f_x[i_ic] - vec_f_x[i_a]))
		fb = float(abs(vec_f_x[i_ic] - vec_f_x[i_b]))
		if max(fa,fb)==0:
			return -999
		return min(fa,fb)/max(fa,fb)

def min_max(vec_f_x, i_nx):
	# finds the minimum and maximum values of vector 'x' (OFFSET = 1).
	# nx: size of 'x'. nx = (long)x[0];

	f_min = vec_f_x[0]
	f_max = vec_f_x[0]

	for i in range(1, i_nx):
		if(vec_f_x[i] < f_min):
			f_min = vec_f_x[i]
		if(vec_f_x[i] > f_max):
			f_max = vec_f_x[i]

	return f_min, f_max

def reverse_minmax(x):
	'''
	transforms vector 'x' in such a way that now the minimum is the
	maximum and conversely (OFFSET = 1).   
	nx: size of 'x'.
	'''
	return x[::-1]

def reverse_minmax_2(x):
	# /*!  <PRE>
	#    transforms vector 'x' in such a way that now the minimum is the
	#    maximum and conversely (OFFSET = 1).
	# Transformation:
	#     x'[i] = -x[i]
	  
	#   nx: size of 'x'.

	#   THE FUNCTION IS DESTRUCTIVE ON (x). IT OVERWRITES (x).
	#   NO MEMORY IS ALLOCATED WITHIN THE FUNCTION.

	return -1*x

def center_max(win, l_WindowLength):

	'''	Determine if the central value of the vector is a local maximum.
	(The function must be monotonously growing before him and the same
	but decreasing after him).
	Input:
	win: vector with the time series (OFFSET = 1).
	l_WindowLength: # of 'win' points.
	Output:
	center = 0 there is no minimum
	1 there is a local maximum
	'''

	lw2  = l_WindowLength//2 #/* greater operator half-length */
	lwt2 = lw2 + 1          #/* greater operator half-length */

	# / * search for maximums. left branch * /
	ok=1
	for i in range (lw2-1, -1, -1): ##i-- in for loop in c code
		if( win[i] > win[i+1] ):  #* non-strict checking inside window */ 
			ok = 0
		if not ok:
			break

	if(ok == 1):
		for i in range(lwt2, l_WindowLength): #/ * search for maximums. right branch * /
			if(win[i] > win[i-1]): #/* non-strict checking inside window */ 
				ok = 0
			if not ok: #ok==0
				break

	if(ok == 1):
		if((win[0] >= win[lwt2]) or (win[l_WindowLength-1] >= win[lwt2])):
			ok=0

	return ok 

def center_min(win, l_WindowLength):
	# / *! <PRE>
	# Determine if the central value of the vector is a local minimum.
	# (The function must be monotonously decreasing before him and the same
	# but growing after him).
	# Input:
	# win: vector with the time series (OFFSET = 1).
	# l_WindowLength: # of 'win' points.
	# Output:
	# center = 0 there is no minimum
	# 1 there is a local minimum

	lw2 = l_WindowLength//2      #/ * lower operator half-length * /
	lwt2 = lw2 + 1        #/ * lower operator half-length * /

	#/ * minimum search. left branch * /
	ok=1
	for i in range(lw2, -1, -1):
		if(win[i] < win[i+1]): #/* non-strict checking inside window */
			ok = 0
		if not ok:
			break

	if(ok == 1):
		for i in range(lwt2+1, l_WindowLength+1): #/ * search for maximums. right branch * /
			if(win[i] < win[i-1]):     #/* non-strict checking inside window */
				ok = 0
			if not ok:
				break

	if(ok == 1):
		if((win[0] <= win[lwt2]) or (win[l_WindowLength] <= win[lwt2])):
			ok=0

	return ok

def centers_list(x, nx, l_ExplorationWindowlLength, i_Minmax):
	'''	
	/ *! <PRE>
	Build the list of points in a series that are maximum
	or relative minimums within an observation window.
	These points are the candidates for anomaly centers.

	NOTE: In order to qualify as an extremum, the function must be monotonote at both sides
	of the candidate (peak) point and also The candidate peak cannot be the next point
	in the signal w.r.t the previous extreme found.
	That is, there must be more than one sampling interval between peaks.

	Input:
	x: vector with the time series (OFFSET = 1).
	nx: size of 'x'.
	l_ExplorationWindowlLength: observation window length (odd)
	i_Minmax = -1 a minimum is searched
	1 maximum is sought

	Output:
	centers_list: vector with the list of found centers
	</PRE> * /
	'''
	l_CurrentNumberOfExtrema=0
	centers=[]
	win=np.zeros((l_ExplorationWindowlLength,))
	lw2  = l_ExplorationWindowlLength//2 #lower operator half length
	npos = nx - l_ExplorationWindowlLength + 1 #number of windows to examine

	for i in range(0, npos):
		l_CandidateExtremumLocation = i + lw2
		for j in range (0, l_ExplorationWindowlLength):   #/* define a window */
			win[j] = x[i + j]

		# if( i_Minmax > 0): #maxs
		j = center_max(win, l_ExplorationWindowlLength)

		if( j == 1 ):
			if(l_CurrentNumberOfExtrema == 0):
				centers.append(l_CandidateExtremumLocation)
				l_CurrentNumberOfExtrema = l_CurrentNumberOfExtrema +1
			else:
				if( (l_CandidateExtremumLocation - centers[l_CurrentNumberOfExtrema-1]) > 1 ):
				 # {/* in order to add an extremum, the function must be monotonote at both sides
					 #  of the candidate point and also it can not be the next point in the signal w.r.t
					 #  the previous extrema found. That is, it must be more than one sampling interval further.
				 #  */
					centers.append(l_CandidateExtremumLocation)
					l_CurrentNumberOfExtrema = l_CurrentNumberOfExtrema +1
	return centers 

def line_segment(x, nx, a, b):
	'''
	Replaces the x values in the interval [a,b] by the
	   linear function in [a, b].
		OFFSET = 1.
	'''

	if(a == b):
		return x

	m = (float(x[b]) - float(x[a]))/(b - a)

	for i in range(a, b+1):
		x[i] = float(m * (i - a) + x[a])
		
	return x

def area_triang(f_x1, f_y1, f_x2, f_y2, f_x3, f_y3):
	'''
	 Calculate the area of a triangle given the coordinates. of the 3
		 vertices
		 Input:
		 x1, y1, etc ... Coordinates of the three vertices of a triangle

		 Output:
		 area: triangle area
	'''
	a = math.sqrt((f_x1-f_x2) * (f_x1-f_x2) + (f_y1-f_y2) * (f_y1-f_y2))
	b = math.sqrt((f_x2-f_x3) * (f_x2-f_x3) + (f_y2-f_y3) * (f_y2-f_y3))
	c = math.sqrt((f_x1-f_x3) * (f_x1-f_x3) + (f_y1-f_y3) * (f_y1-f_y3))
	s = (a + b + c)/2
	s = s * (s-a) * (s-b) * (s-c)

	if(abs(s) <= 1e-10):
		area = 0
	else:
		area = math.sqrt( s )

	return area

def clip(vec_f_x, i_nx, i_i, xi, yi):
	'''
	#################
	#check points.... i1 and i2
	#################
	clips the coordinate of a point from a time serie wrt
	the initial and final points
	Input
	x: series (OFFSET = 1).
	nx: size of 'x'.
	i: x-coordinate
	Output
	xi, yi: x, y value associated to the given point
	'''
	i1 = 0
	i2 = i_nx

	if(i_i < i1):
		xi=float(i1)
		return xi,yi
	if(i_i >= i2):
		xi = float(i2)
		return xi,yi

	xi = float(i_i)
	yi = vec_f_x[int(xi)]

	return xi,yi

def even(i_n):
	#return 1 if even 0 is odd
	if( (i_n % 2) == 0):
		return 1
	else:
		return 0

def area_pico(vec_f_x, i_nx, i_ic, f_threshold, i_ii, i_id):
	'''
	 From the position of a supplied center (ic) of a
		 given series (x) produces the sequence of areas under the peak of
		 MAXIMUM up to background level
		 Input:
		 x: vector with the time series (OFFSET = 1).
	  nx: size of 'x'.
		 ic: point corresponding to the center of the anomaly to study
		 threshold: threshold for the ratio of two successive areas "going down"
				 from the position of the anomaly that is in 'ic'.
				 When the ratio between two succesive areas falls under the threshold
				 the expansion stops.

		 Output:
		 area: area of the peak corresponding to 'ic'
		 ii, id: abscissa that define the beginning and end of the peak
	'''
	x1,y1,x2,y2,x3,y3 = None, None, None, None, None, None
	x1, y1 = clip( vec_f_x, i_nx, i_ic-1, x1, y1)
	x2, y2 = clip( vec_f_x, i_nx, i_ic  , x2, y2)
	x3, y3 = clip( vec_f_x, i_nx, i_ic+1, x3, y3)

	aold = area_triang(x1, y1, x2, y2, x3, y3)
	iizq = i_ic - 1
	ider = i_ic + 1
	stopi = 0
	stopd = 0
	stopArea = 0
	l_MaximumNumberOfTrials = 0
	j = 1

	while True:
		j = j+1
		if( even(j) and (stopi == 0)):    #move fwd by one step on the left
			x2, y2 = clip(vec_f_x, i_nx,   iizq, x2, y2)
			iizq = iizq-1
			x1, y1 = clip(vec_f_x, i_nx, iizq, x1, y1) ###### --iizq
			x3, y3 = clip(vec_f_x, i_nx,   ider, x3, y3)
			
			i = int(x2)
			iiold = i
			k = int(x3)
			idold = k
			if(y1 > y2):    #/* slope reversal */
				# /* go back */
				iizq = iizq+ 1
				iiold = iizq
				stopi = 1 #   /* congelar el punto en lugar de retornar */
			
		if(even(j) == 0 and stopd == 0):    #/* avance de 1 punto por la derecha */
			x1, y1 = clip(vec_f_x, i_nx,   iizq, x1, y1)
			x2, y2 = clip(vec_f_x, i_nx,   ider, x2, y2)
			ider=ider+1
			x3, y3 = clip(vec_f_x, i_nx,   ider, x3, y3) ###### ++ider
			i = int(x1)
			iiold = i
			k = int(x2)
			idold = k
			if(y3 > y2): #    /* slope reversal */
				ider = ider-1
				idold = ider
				stopd = 1
		delta_a = area_triang(x1, y1, x2, y2, x3, y3)

		i_ii = iiold
		i_id = idold
		anew = aold + delta_a

		if (aold != 0. ):
			ratio = anew/ aold
		else:
			ratio = float(1.0e18)

		if(ratio < f_threshold):
		  stopArea = 1

		aold = anew
		l_MaximumNumberOfTrials = l_MaximumNumberOfTrials +1
		if ( (l_MaximumNumberOfTrials <= i_nx) and (not stopArea) and (stopi == 1 and stopd == 1) ):
			break

	return aold, i_ii, i_id





