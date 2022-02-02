# PioneerMATLAB
Respository of MATLAB functions for use with BORL Pioneer system

SL: I wrote help comments for these functions, so you can just do the usual 'help find_intensity' and it will tell you the format of the inputs.

List of functions:

find_intensity - gives an array with the sum of the histogram at each pixel
find_m1 - gives an array with the computed first moment (mean) at each pixel, need to specify a threshold (e.g. 0.05 for 5% of the peak)
nan_mean - computes mean of a vector ignoring NaN entries (mean would just return NaN as soon as there is one NaN entry)
nan_sum - computes sum of a vector ignoring NaN entries
plot_timing - plots timing response of a section of an array (e.g. can be column, row, sub-section or whole array)
surf_intensity - does a surf of the intensity of a timing_response input (R x C cell)
surf_m1 - does a surf of the m1 of a timing_response input (R x C cell)
