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

waveList = [ 725 838 ];   % Example [689 725 802];

sources =  [7];

path_prj = fullfile(Path_base, proj_name);
hdf5_file = fullfile(path_prj, strcat(Prefix,'-', proj_name,'.hdf5'));
% hdf5_file = [Path_base '/' proj_name '/laser-20220228_PETG_homo.hdf5'];
info = h5info(hdf5_file, '/Data');

if exist(fullfile(path_prj, name_fld_corrected), 'dir')==0
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
    
    for col_i = 1:cols

        for row_i = 1:rows
            
            if isnan(add_corr(row_i,col_i)) || isempty(timing_response(:,row_i,col_i)) ||...
                   ((row_i==17&&col_i==24)) || (col_i == 32) && (row_i == 2)  % changed by jingjing 2021.10.11
                timing_response(:, row_i, col_i) = 0; % added by jingjing 2021.05.11
            else
                timing_response(:, row_i, col_i) = correct_histogram(timing_response(:, row_i, col_i), add_corr(row_i, col_i), col_cal(col_i));
             
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