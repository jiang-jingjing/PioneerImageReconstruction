function my_mean = nan_mean(vector)
%%% calculates mean of vector where all NaNs are excluded

vector(isnan(vector)) = [];
my_mean = mean(vector);

end