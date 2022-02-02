% Script to plot intensity from timing response cell for arbitrary sized array
% Function of the form intensity = find_intensity(timing_response, rows, cols)
% timing_response - cell where each element is the timing response of a
% pixel.
% rows - range of rows to include in the intensity calculation, e.g. 7:19
% cols - range of cols to include in the intensity calculation, e.g. 12:24
%
% This script subtracts dark counts!!

function intensity = surf_intensity(timing_response, rows, cols)

intensity = nan(length(rows),length(cols));

for i_row = 1:length(rows)
    for i_col = 1:length(cols)
        row = rows(i_row);
        col = cols(i_col);
        if isnan(timing_response{row,col})
        else
            x = timing_response{row,col};
            x = abs(x - median(x));
            intensity(i_row,i_col) = sum(x);
        end
    end
end

intensity(2,32) = 0;

figure;
surf(intensity);
xlabel('Column');
ylabel('Row');

end