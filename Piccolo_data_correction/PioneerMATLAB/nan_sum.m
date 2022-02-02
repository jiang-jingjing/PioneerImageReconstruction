function my_sum = nan_sum(vector)
%%% calculates mean of vector where all NaNs are excluded

vector(isnan(vector)) = [];
my_sum = sum(vector);

end