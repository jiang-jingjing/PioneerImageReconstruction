% Script to plot intensity from timing response cell for arbitrary sized array
% Function of the form intensity = find_intensity(timing_response, rows, cols)
% timing_response - cell where each element is the timing response of a
% pixel.
% rows - range of rows to include in the intensity calculation, e.g. 7:19
% cols - range of cols to include in the intensity calculation, e.g. 12:24
%
% This script subtracts dark counts!!
% modified April, 2021 by jingjing jiang
function intensity = find_intensity(timing_response, rows, cols)

if  iscell(timing_response)
intensity = nan(length(rows),length(cols));

for i_row = 1:length(rows)
    for i_col = 1:length(cols)
        row = rows(i_row);
        col = cols(i_col);
        if isnan(timing_response{row,col})
        else
            x = timing_response{row,col};
%             x = abs(x - median(x));
            intensity(i_row,i_col) = sum(x);
        end
    end
end

if length(cols) > 31
    intensity(2,32) = 0;
end
else
intensity = squeeze(sum(timing_response,1));
end
end