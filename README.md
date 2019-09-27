# shtfinder

module shtfinder for search in Globus-M plasma pulses
made for calculation the euclidian and DTW distance of the 
plasma discharge diagnostics time series and consecutive sorting

Created by V. Solokha
Last updated    2019-09-07: Function created
                2019-09-20: WINE used for linux usage of test.exe 
                2019-09-27: Improved input reading routine, merged with convert.py

==========
MODULES:
==========

pandas, numpy, matplotlib, scipy, fastdtw, tqdm

==========
USAGE:
==========
python3 shtfinder.py baseshot start_sequence end_sequence diagnostic_idx weights shtpath

Example: python3 36612 36610 36614 6,7 1.0,0.0 "./sht"

==========
INPUTS:
==========

baseshot (int)         -- the number of the base pulse
start_sequence (int)   -- starting point of the range
end_sequence (int)     -- finishing point of the range 
diagnostic_idx (array) -- indices of the comparable diagnostics 
weights (array)        -- weights of the diagnostics
shtpath (str)          -- path to the sht files folder
