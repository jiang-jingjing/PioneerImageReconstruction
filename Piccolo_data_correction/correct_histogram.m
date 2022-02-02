function [histogram] = correct_histogram(x, add_corr, col_cal) %#codegen
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%%%% Load signal and correct for dnl, added 12.11.18 SL
sig_interp = 100;

x = x - median(x);  %% Subtracts any dark counts

x_counts = sum(x);

%%% Need to apply any correction of DNL before
%%% interpolation etc!! %%%

x = flipud(x);
% LM: interp1 seems much faster than interp, though I am not sure about the
% accuracy tradeoff
y = interp1(0:(length(x)-1), double(x),(0:(length(x)*sig_interp)-1)/sig_interp, 'spline');   %% interpolate x at rate sig_interp
    %interp(double(x), sig_interp);
l_start = length(y);        %% length before scaling

%l_end = length(y)*col_cal(col_i);   %% length after scaling

%%%% This part of the code just makes y the correct length,
%%%% this needs to account for the slicing of the TPSF

l_end = (256*sig_interp)/col_cal;

if l_end > length(y)
    sig_diff = round(l_end - length(y));
    y = [y(1:end); zeros(sig_diff,1)];
end

if l_end < length(y)
    sig_diff = round(length(y) - l_end);
    y = y(1:(end - sig_diff));
end

%%%%% Add shift for additive part of the chip calibration

bins_shift = (add_corr*sig_interp)/col_cal;
bins_shift = round(bins_shift);

%try
    y = circshift(y,-bins_shift);
%catch
%    temp = 1;
%end

%%%%% Interpolate to 256 bins to correct for LSB

x1 = 1:(length(y));
x2_start = 1; % starts at 1 because of the flip.
x2_end = (length(y)-(length(y)/256));
x2 = x2_start:((x2_end-x2_start)/255):x2_end;
y1 = interp1(x1, y, x2); % spline is unneccessay when downsampling LM
y1 = round(y1)'; % transpose for convenience
TXX = y1(:) < 0;
y1(TXX) = 0;

%%%%% Finally correct number of photons to initial value so
%%%%% that intensity is correct

y_counts = sum(y1); %% Total counts after correction

for i_y = 1:length(y1)
    y1(i_y) = (y1(i_y)/y_counts)*x_counts;
end

%%%%% And load back into timing response
histogram = flipud(y1);
end

