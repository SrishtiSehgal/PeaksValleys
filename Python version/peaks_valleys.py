import numpy as np
import peaks_valleys_lib as lib
import matplotlib.pyplot as plt

def load_signal(str_file, column=(0,), i_skiprow=1, delimiter=','):
	return np.loadtxt(str_file, skiprows=i_skiprow, delimiter=',', usecols=column)

def make_minima_or_maxima_database(vf1_Signal, l_SignalLength, f_xOfTheFirstPoint,
	dx, i_Minmax, l_WindowLength, c_ExtremumFindingCriterion, f_Threshold):
	'''
		From a signal float builds a database in the form of a
		RAF file
		Input:
		nparr_signals: senyal (float) to process. As always, in [0] is the number
				 Of elements.
		int_SignalLength: No. of signals of the signal
		int_xOfTheFirstPoint: x-value corresp. to the first element of the signal.
		int_dx: sampling interval.
		int_Minmax: Indicator of the type of object to be detected:
				-1 minimum
				 1 max
				 0 both
		int_WindowLength: Length of the observation window for the detection of the maximums/minimums ODD number
		char_ExtremumFindingCriterion: Approach used for finding the atart and end points of a
								maximum / minimum: {"area", "derivative"}.
								area: The maximum is extended one point at a time away from the
									  central location and the area underneath is compared
									  with the area prior to the expansion.
									  When the ratio between two succesive areas falls under
									  the threshold the expansion stops.
								derivative: The (finite difference) derivative at both sides of the
											extremum is computed. If a given derivative falls
											under the threshold or reverses sign, the expansion in the given
											direction stops (it may continue in the other direction).

		f_Threshold: If (s_ExtremumFindingCriterion) is 'area', threshold for the ratio of areas along the end that is scanned.
					 If (s_ExtremumFindingCriterion) is 'derivative', threshold for the directional derivative. It the sign reverses or
													 if slope falls under this value, the expansion along that direction stops.
					 Required to determine the lateral extension of the maximum / minimum
		(Ex: 0.001).
		str_DatabaseMatrixOutputFile: Name of the database matrix output file (NO EXTENSION).
									IF NULL NO OUTPUT FILE IS CREATED
		str_FileType: {"txt", "csv"}.
		txt: text file.
		csv: comma delimited file with a header composed of the computed field values ​​(see below).

		database: matrix with the database (offset = 1).
							  database [0] [0] contains the number of records
							  database [0] [1] contains the number of fields
		Ouput:
		err = 0 all OK.
		= -1 The database could not be created.
		= -2 'i_Minmax' wrong.
		= -3 'l_WindowLength' out of range.
		= -4 negative 'threshold_deriv'.
		= -5 No end found.
		= -6 wrong criterion for extremum finding (not "area", now "derivative").

		Format of the database created:
		db [0] [0]: number of extreme
		db [0] [1]: number of fields (16)
		Each record is a vector (db) of 16 elements with the following variables:
		db [1] = (float) index;
		db [2] = (float) position;
		db [3] = value;
		db [4] = amplitude;
		db [5] = (float) start;
		db [6] = (float) end;
		db [7] = (float) width;
		db [8] = xPosition; The same as [2], [5], [6] and [7] but in x-coordinates.
		db [9] = xStart;
		db [10] = xEnd;
		db [11] = xWidth;

		db [12] = x_symmetry;
		db [13] = y_symmetry;
		db [14] = sigma_left;
		db [15] = sigma_right;
		db [16] = (float) i_Minmax; 1 = maximum, -1 = minimum

		File "csv":
		1) header with variable names
					index
		position
		value
		amplitude
		start
		end
		width
					xPosition
					xStart
					xEnd
					xWidth
		x_symmetry
		y_symmetry
					sigma_left
					sigma_right
		i_Minmax
		2) data matrix with the info about the extreme (min / max) found.
		*vf1_LinearizedBackground: signal where the extrema (maxima, minima) have been replaced by linear 
								  segments connecting their bases (like cutted with a scissors).
								 
		*vf1_Residual            : Residual signal when (*vf1_LinearizedBackground) is removed from
								  the original.

		Call:
		int make_minima_or_maxima_database
							( vf1_Signal, l_SignalLength
							, f_xOfTheFirstPoint, f_dx
							, i_Minmax, l_WindowLength
							, c_ExtremumFindingCriterion
							, f_Threshold
							  , s_DatabaseMatrixOutputFile
							, s_FileType
							, &mf1_OutputDatabase
								 , &l_NumberOfAttributesOfTheDatabase
							, &l_NumberOfExtremaFound
							, &vf1_LinearizedBackground
							, &vf1_Residual
							);
	'''
	#############################
	#initialize

	err = 0
	i = 0
	clist=None
	db=None
	start=0
	end=0

	#############################
	#error checking
	if( (vf1_Signal is None) or (l_SignalLength == 0)):
		return -1
	if( (i_Minmax != -1) and (i_Minmax!= 1)):
		return -2
	if( (l_WindowLength < 0) or (l_WindowLength > l_SignalLength)):
		return -3
	if(f_Threshold < 0):
		return -4

	#############################
	# /* specifying the appropriate extremum finding function */
	if(c_ExtremumFindingCriterion == 'a'): #/* "area"*/
		extremumFindingCriterionFunction = lib.area_pico
	else:
		if(c_ExtremumFindingCriterion == 'd'): # /* "derivative" */
			extremumFindingCriterionFunction = lib.peak_amplitude
		else:
			return -6

	#############################
	# Finds the global min and max of the signal
	# f_min, f_max = min_max(vf1_Signal, len(vf1_Signal), f_min, f_max)
	
	# number of variables to calculate for each anomaly		
	# nvar = numberOfFields

	# /* construction of (vf1_LinearizedBackground) */
	vf1_LinearizedBackground = vf1_Signal
	vf_WorkingSignal = np.copy(vf1_Signal)

	#if the user requires to find minimums, flip signal
	if( i_Minmax == -1):
		vf_WorkingSignal = lib.reverse_minmax_2(vf_WorkingSignal)

	# /* DETECTION OF THE EXTREMA POINTS (maxima, minima) */
	clist = lib.centers_list(vf_WorkingSignal, len(vf_WorkingSignal), l_WindowLength, 1)
	if(clist is None):
	   return -5

	mf1_OutputDataBase = None
	user_break, error_peak=False, False
	l_NumberOfAttributesOfTheDatabase = 16

	for i in range(0, len(clist)):
		db = np.zeros((1, l_NumberOfAttributesOfTheDatabase))

		position      = clist[i] #positioning starts at 0, array index-like
		if i_Minmax<0:
			value = -(vf_WorkingSignal[clist[i]])
		else:
			value = vf_WorkingSignal[clist[i]]

		amplitude, start, end = extremumFindingCriterionFunction(vf_WorkingSignal, l_SignalLength, clist[i], f_Threshold, start, end)

		width		  = (end-start)
		xPosition     = (position-1)* dx + f_xOfTheFirstPoint
		xStart        = (start   -1)* dx + f_xOfTheFirstPoint
		xEnd          = (end     -1)* dx + f_xOfTheFirstPoint
		xWidth        = (width     )* dx

		if start==end:
			error_peak = True
			user_break=False
			print("{} has errors. Seems that immediate values around it are the same as the peak.".format(position))
			if input('Stop processing? Hit the "B" key.')=='B':
				user_break=True
				return -7
		
		if not user_break:
			x_symmetry    = lib.simetria_pico(position, start, end)
			y_symmetry    = lib.simetria_valores(vf_WorkingSignal, position, start, end)
			sigma_left    = float(lib.sigma2(  value,  xPosition,  xStart, vf_WorkingSignal[start]))
			sigma_right   = float(lib.sigma2(  value,  xPosition,  xEnd  , vf_WorkingSignal[end  ]))
	
			db[ 0,0] = i
			db[ 0,1] = float(position)
			db[ 0,2] = value
			db[ 0,3] = amplitude
			db[ 0,4] = float(start)
			db[ 0,5] = float(end)
			db[ 0,6] = float(width)
				   
			db[ 0,7] = xPosition
			db[ 0,8] = xStart
			db[ 0,9] = xEnd
			db[ 0,10] = xWidth
	
			db[ 0,11] = x_symmetry
			db[ 0,12] = y_symmetry
			db[ 0,13] = sigma_left
			db[ 0,14] = sigma_right
			db[ 0,15] = float(i_Minmax)
	
			# /* computation of the i-th linear background segment */
			vf1_LinearizedBackground = lib.line_segment(vf1_LinearizedBackground, l_SignalLength, start,  end)
	
			if mf1_OutputDataBase is None:
				mf1_OutputDataBase = db
			else:
				mf1_OutputDataBase = np.vstack((mf1_OutputDataBase, db))
	
	# /* residual */
	vf1_Residual = vf_WorkingSignal - vf1_LinearizedBackground

	# /* [1]: position [3]: start [4]: end */ ################################
#	db[0,1] = f_xOfTheFirstPoint + (db[0,1]-1) * dx#; /* correction of the peaks position, start and end points according to the */
#	db[0,3] = f_xOfTheFirstPoint + (db[0,3]-1) * dx#; /* sampling interval and the x-value of the first point. */
#	db[0,4] = f_xOfTheFirstPoint + (db[0,4]-1) * dx

	if clist is not None:
		# /* export the resulting matrix 
		np.savetxt('Peaks and valleys output.csv', 
				mf1_OutputDataBase, 
				header= "index,position,value,amplitude,start,end,width,xPosition,xStart,xEnd,xWidth,x_symmetry,y_symmetry,sigma_left,sigma_right,i_Minmax", 
				delimiter=',', 
				comments='')
		
		np.savetxt('vf1_LinearizedBackground_output.csv', 
				   vf1_LinearizedBackground, 
				   delimiter=',', 
				   comments='')
		
		np.savetxt('vf1_Residual_output.csv', 
				   vf1_Residual, 
				   delimiter=',', 
				   comments='')
		t = range(len(vf_WorkingSignal))
		plt.plot(t, vf1_Signal, 'r', t, vf1_LinearizedBackground, 'b', t, vf1_Residual, 'g')
		plt.show()
	else:
		clist = None
		vf_WorkingSignal = None
	if error_peak:
		return -7
	else:
		return(err)


if __name__ == '__main__':
	#	vf1_Signal = load_signal(input('filepath: '), 
	#		tuple(input('comma separated list of columns to read e.g. 0, or 1,2,3: ')), 
	#		int(input('rows to skip: ')))
	
	vf1_Signal = load_signal('xxx-original-signal-1.csv', (1,), 1) 

	int_window_length = 11
	i_Minmax = 1
	c_ExtremumFindingCriterion = 'd'

	f_Threshold = 0.0001
	int_xCoordinateOfTheFirstPoint = 0
	int_sampling_interval = 1

#	int_window_length = int(input('length of the operator for finding an extremum (ODD number): '))
#	i_Minmax = -1 if input('Type of extrema {minima, maxima}: ')=='minima' else 1
#	c_ExtremumFindingCriterion = str('Type of extrema criterion {area=a, derivative=d: ')
#
#	print('Stopping criterion for finding the limits of an extremum at the left and right sides.')
#	print('If \"area\", a value in [0,1] as an area ratio. Typically ~0.98 works well.')
#	print('If \"derivative\", a threshold for the derivative value.')
#	print('')
#
#	f_Threshold = float(input('Enter proper threshold value: '))
#	int_xCoordinateOfTheFirstPoint = int(input('x coordinate of the first point. Default 0: '))
#	int_sampling_interval = int(input('sampling interval. Default is 1: '))
	
	print('Error code information:')
	print(' 0 all OK.')
	print('-1 The database could not be created.')
	print('-2 i_Minmax wrong.')
	print('-3 l_WindowLength out of range.')
	print('-4 negative threshold_deriv.')
	print('-5 No end found.')
	print('-6 wrong criterion for extremum finding (not "area", now "derivative").')
	print('-7 data base values for one or more of the peaks were incorrect. Please check log.')
	print('')
	
	int_error_code = make_minima_or_maxima_database(vf1_Signal, vf1_Signal.shape[0], int_xCoordinateOfTheFirstPoint, 
												 int_sampling_interval, i_Minmax, int_window_length, c_ExtremumFindingCriterion, f_Threshold) 

	print('Error code produced: ', int_error_code)