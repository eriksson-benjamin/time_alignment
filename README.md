# time_alignment
Perform time-alignment of TOFu DAQ.

Generate a data set using create_TOF.py with the --save-data argument. Store data sets in '/common/scratch/beriksso/TOFu/data/time_alignment/'.

An example input file for create_TOF.py used for dataset_2 is shown below:

--disable-bgs
--disable-cuts
--disable-plots
--save-data /common/scratch/beriksso/TOFu/data/time_alignment/dataset_2/
--time-range-file input_files/time_ranges.txt
--disable-scratch

calculate_coincidences.py
-------------------------
Calculates coincidences between S1-5 and S2 detectors. Selects the events contirbuting to the gamma peak and uses these to calculate coincidences between S1-5 and the other S1's. Coincidences are saved to file stored in output/.


