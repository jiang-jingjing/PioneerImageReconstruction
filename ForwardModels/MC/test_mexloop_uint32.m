% test mexloop_uint32.m
% step 1: in MATLAB command window, type in: mex mexloop_uint32.c 
%         to compile mexloop_uint32
% step 2: replace the for-loop with HIST = mexloop_uint32(len, detID, bin, Track_atten, N_DET, N_BIN);

% ---- dummy data to test the function ----%
len = 10;
detID = [1 2 2 2 3 1 2 2 2 3]';
bin = [1 2 3 1 4 1 2 3 1 4]';
Track_atten = [1:len];
N_DET = max(detID);
N_BIN = max(bin);
% ---- dummy data to test the function ----%

 
detID = uint32(detID);
bin = uint32(bin);
N_DET = uint32(N_DET);
N_BIN = uint32(N_BIN);

Track_atten = single(Track_atten); % Skip this if Track_atten is single 

HIST = mexloop_uint32(len, detID, bin, Track_atten, N_DET, N_BIN);
