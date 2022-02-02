% Script to plot timing response for a timing response over rows and cols,
% timing response - cell of arbitrary size, where each cell has and N x 1 
% timing response.
% rows - specifies range of rows, e.g. 5:20,
% cols - specifies range of cols, e.g. 5:20

function [timing_plot] = plot_timing(timing_response, rows, cols)

timing_plot = figure;
hold on
for i_row = 1:length(rows)
    for i_col = 1:length(cols)
        row = rows(i_row);
        col = cols(i_col);
        if isnan(timing_response{row,col})
        else
        x = timing_response{row,col};
        plot(x);
        end
    end
end

end


