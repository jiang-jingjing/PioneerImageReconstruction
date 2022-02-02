% Analysis script to plot first moment of timing response arbitrary sized 
% array, m1 = find_m1(timing_response, threshold, rows, cols)
% timing response - cell of abitrary size, where each cell location
% contains a timing response N x 1 bins.
% threshold - defines the threshold at which to exclude bins from the
% calculation, e.g. 0.05, sets a threshold at 5% of the response maximum.
% rows - defines range of rows to comput timing response over, e.g. 1:5
% cols - defines range of colss to comput timing response over, e.g. 1:5

function m1 = surf_m1(timing_response, threshold, rows, cols)

m1 = nan(length(rows),length(cols));

for i_row = 1:length(rows)
    for i_col = 1:length(cols)
        row = rows(i_row);
        col = cols(i_col);
        if isempty(timing_response{row,col}) || (row == 2 && col == 32)
        else
        tr1 = timing_response{row,col};
        [maxC, maxI] = max(tr1);
        A = tr1 > maxC*threshold;
        A = A';
        m1(i_row,i_col) = sum((transpose(tr1).*(1:length(tr1)).*A))/sum(transpose(tr1).*A);
        end
    end
end

figure;
surf(m1)
xlabel('Column');
ylabel('Row');

end