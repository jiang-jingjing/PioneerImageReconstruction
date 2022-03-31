%%%%%%%% Timing data correction script so that we can do wrapping %%%%%%%%%

%%% Modified the script March 2022

%%%% Put the timing responses that you want to correct into the timing_data
%%%% folder, add the name of the files to the source and how many sources
%%%% you want to correct, then run!

clear all;
format compact;

%ParPool creation
numcores = feature('numcores'); %Get number of physical cores
p = gcp('nocreate');

if isempty(p)
    parpool('local', numcores); %create a parpool if doesn't exist
elseif p.NumWorkers < numcores
    delete(p)
    parpool('local', numcores); %re-create the pool if the numcore used is not the maximum possible
end

%addpath('PioneerMATLAB/')


Path_base = ['/media/jiang/WD10T/Data/Projects/PioneerDebug2022/Data/'];

proj_name = 'ICGlipid_20220311';

name_fld_corrected = 'timing_data_corrected';
Prefix = 'laser';

waveList = [ 780 830 ];   % Example [689 725 802];

sources =  [1:11];

path_prj = fullfile(Path_base, proj_name);
hdf5_file = fullfile(path_prj, strcat(Prefix,'-', proj_name,'.hdf5'));
% hdf5_file = [Path_base '/' proj_name '/laser-20220228_PETG_homo.hdf5'];
info = h5info(hdf5_file, '/Data');

if exist(fullfile(path_prj, name_fld_corrected))==0
    mkdir(fullfile(path_prj, name_fld_corrected))
else
    rmdir(fullfile(path_prj, name_fld_corrected), 's') %remove it first if already exist
    mkdir(fullfile(path_prj, name_fld_corrected))
end

tic
% find all datasets to process
ids = zeros(1, length(info.Datasets));
file_pos_ls = strings(1,length(info.Datasets));

parfor index = 1:length(info.Datasets)
    if ismember(h5readatt(hdf5_file, ['/Data/' num2str(index)], 'Wavelength [nm]'), waveList) && ismember(h5readatt(hdf5_file, ['/Data/' num2str(index)], 'Channel'), sources)
        
        ids(index) = index;
        
        datapath = ['/Data/' num2str(index)];

        target_pos = h5readatt(hdf5_file, datapath, 'Targetted Head position {x, y, z, dist}');
        file_pos_ls(index) = ['target_' num2str(target_pos(1)) '_' ...
        num2str(target_pos(2)) '_' ...
        num2str(target_pos(3)) '_' ...
        num2str(target_pos(4)) ];

    end
end
fprintf('Time elapsed for file reading: %f sec\n',toc);

ids = nonzeros(ids);
file_pos_ls = unique(file_pos_ls);

%Create target_pos folder outside parfor
for file_pos = file_pos_ls
    if ~strcmp(file_pos,"") && exist(fullfile(path_prj, name_fld_corrected, file_pos), 'file')==0
        mkdir(fullfile(path_prj, name_fld_corrected, file_pos))
    end
end

%% check powermeter

%Var names needed in read_pw
file = hdf5_file;
fldr = convertStringsToChars(path_prj);

read_pw

%%
tic
%%%% Load column calibration
temp = load('camera_correcting/col_calibration.mat');  
col_cal = temp.col_calibration;

%%%%%% Load additive calibration
temp = load('camera_correcting/add_corr.mat');
add_corr = temp.add_corr;

rows = 32;
cols = 32;

parfor i_id = 1:length(ids)
    %%% Load data
    datapath = ['/Data/' num2str(ids(i_id))];
     
    timing_response = double(h5read(hdf5_file, datapath));
    wav = num2str(h5readatt(hdf5_file, datapath, 'Wavelength [nm]'));
    source_no = h5readatt(hdf5_file, datapath, 'Channel');
    rep_no = h5readatt(hdf5_file, datapath, 'Run');
    target_pos = h5readatt(hdf5_file, datapath, 'Targetted Head position {x, y, z, dist}');
    
    file_pos = ['target_' num2str(target_pos(1)) '_' ...
        num2str(target_pos(2)) '_' ...
        num2str(target_pos(3)) '_' ...
        num2str(target_pos(4)) ];

    %sig_interp = 100;
    %timing_response_1 = timing_response;
    %timing_data_tmp = zeros(size(timing_response)); % added by jingjing 2021.05.11
    
    %Generate static variables once
    shared_data = struct;
    shared_data.sig_interp = 100;
    shared_data.t1 = 0:(length(timing_response)-1);
    shared_data.t3 = (0:(length(timing_response)*shared_data.sig_interp)-1)/shared_data.sig_interp;
    
    for col_i = 1:cols
                
        shared_data.l_end = (256*shared_data.sig_interp)/col_cal(col_i);

        for row_i = 1:rows
            
            if isnan(add_corr(row_i,col_i)) || isempty(timing_response(:,row_i,col_i)) ||...
                   ((row_i==17&&col_i==24)) || (col_i == 32) && (row_i == 2)  % changed by jingjing 2021.10.11
                timing_response(:, row_i, col_i) = 0; % added by jingjing 2021.05.11
            else
                timing_response(:, row_i, col_i) = correct_histogram_par(timing_response(:, row_i, col_i), add_corr(row_i, col_i), col_cal(col_i), shared_data);
             
            end
            
        end
    end
    
    foo_name = sprintf('timing_response_%s_%s_%s_%s.mat', Prefix, wav, num2str(source_no), num2str(rep_no));
    save_name = fullfile(path_prj, name_fld_corrected, file_pos, foo_name);   
    
    %timing_response(isnan(timing_response)) = 0;
    
    save_correct_hist(save_name, timing_response);
    
end

fprintf('Time elapsed for computation: %f sec\n',toc);


function save_correct_hist(file_name, timing_response)
        save(file_name, 'timing_response');
end


function [histogram] = correct_histogram_par(x, add_corr, col_cal, data_struct)
%%%Modified version of correct_histogram.m that avoid some redudant computation
%%%%% Load signal and correct for dnl, added 12.11.18 SL
sig_interp = data_struct.sig_interp;

x = x - median(x);  %% Subtracts any dark counts

x_counts = sum(x);

%%% Need to apply any correction of DNL before
%%% interpolation etc!! %%%

x = flipud(x);
% LM: interp1 seems much faster than interp, though I am not sure about the
% accuracy tradeoff

y = interp1(data_struct.t1, double(x), data_struct.t3, 'spline');   %% interpolate x at rate sig_interp

%l_start = length(y);        %% length before scaling

%l_end = length(y)*col_cal(col_i);   %% length after scaling

%%%% This part of the code just makes y the correct length,
%%%% this needs to account for the slicing of the TPSF

l_end = data_struct.l_end;

if l_end < length(y)
    sig_diff = round(length(y) - l_end);
    y = y(1:(end - sig_diff));
    
elseif l_end > length(y)
    sig_diff = round(l_end - length(y));
    y = [y(1:end); zeros(sig_diff,1)];
    
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
