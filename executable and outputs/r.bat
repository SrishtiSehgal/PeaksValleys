::
:: anomaly-for-iar           (jjv)
:: USE:\nanomaly-for-iar <fileWithTwoSignals> <window length (odd)> <-1: minima, 1: maxima>
::                       <ss_ExtremumFindingCriterion> <threshold> 
::                       <xCoordinateOfTheFirstPoint> <sampling interval> 
:: 
:: where:
:: fileWithTwoSignals         : Two-column file with the signals (column vectors): load signal and the target.
:: window length              : length of the operator for finding an extremum (ODD number).
:: minimaORmaxima             : Type of extrema {minima, maxima}.
:: ss_ExtremumFindingCriterion: {\"area\", \"derivative\"}.
:: threshold                  : Stoping criterion for finding the limits of an extremum
::                              at the left and right sides.
::                              If \"area\", a value in [0,1] as an area ratio.
::                                  Tipically ~0.98 works well.
::                              If \"derivative\", a threshold for the derivative value.
::                                  It is hard to recol_mmend values in this case, because of
::                                  the scale of the input signal.n\
:: BaseNameForOutputFiles     : Prefix for naming the output files (NO EXTENSION).
:: xCoordinateOfTheFirstPoint : if not given it defaults to 1.
:: sampling interval          : if not given it defaults to 1.
:: 
::e:\valdes\applications\anomaly-for-iar\proj\vc2010\anomaly-for-iar\x64\Release\anomaly-for-iar-64-r.exe signals.csv 11 maxima area 0.98  xxx
anomaly-for-iar-64-r.exe xxx-original-signal-1.csv 11 maxima area 0.98 xxx
