# shtfinder

module shtfinder for search in Globus-M plasma pulses
made for calculation the euclidian and DTW distance of the 
plasma discharge diagnostics time series and consecutive sorting

### CHANGES:
* 2019-12-08 Added GUI based on Kivy library
* 2019-12-11 Added pure Python sht file decoder written by N.Zhiltsov
* 2021-11-28 Minor bug fix
### REQS:

pandas, numpy, matplotlib, scipy, fastdtw, tqdm, kivy, shtreaper

wine for linux

### USAGE:
python3 src/gui.py

or

python3 shtfinder.py baseshot start_sequence end_sequence diagnostic_idx weights shtpath

Example: python3 shtfinder.py 36612 36610 36614 6,7 1.0,0.0 "./sht"

### INPUTS:

* baseshot (int)         -- the number of the base pulse
* start_sequence (int)   -- starting point of the range
* end_sequence (int)     -- finishing point of the range 
* diagnostic_idx (array) -- indices of the comparable diagnostics 
* weights (array)        -- weights of the diagnostics
* shtpath (str)          -- path to the sht files folder
