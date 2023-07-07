%%%%%%%% Timing data correction script so that we can do wrapping %%%%%%%%%

%%% Modified the script April 2021

%%%% Put the timing responses that you want to correct into the timing_data
%%%% folder, add the name of the files to the source and how many sources
%%%% you want to correct, then run!
clear all;
format compact;
tic
%temp = load('./dnl_inl_v2/dnl.mat'); % Added by SL 12.11.2018
%dnl = temp.dnl;
addpath('PioneerMATLAB/')
Prefix = 'laser';

waveList = [725];

sources =  [1:11];
% proj_name = 'SPih10r20230503'
FOLDER = ['/media/jiang/WD10T/Data/Projects/PioneerExperiment2023/Data/']
flag_prj = ['20230616_TS_lq5_1'];
fldr = [FOLDER flag_prj];
file = [fldr '/laser-' flag_prj '.hdf5'];
% flag_prj = '/lp_20210615_depth2incl/lq_dist_11_6/';
% fldr = [FOLDER flag_prj];
% file = [fldr '/laser-lq_dist_11_6.hdf5'];

info = h5info(file, '/Data');
if exist([fldr '/timing_data_corrected'])==0
    mkdir([fldr '/timing_data_corrected'])
end
 
% find all datasets to process
ids = [];
for index = 1:length(info.Datasets)
    if ismember(h5readatt(file, ['/Data/' num2str(index)],...
            'Wavelength [nm]'), waveList) &&...
            ismember(h5readatt(file, ['/Data/' num2str(index)], 'Channel'), ...
            sources)
        ids = [ids index];
    end
end

%% check powermeter
read_pw

%%
%%%% Load column calibration
temp = load('camera_correcting/col_calibration.mat');  
col_cal = temp.col_calibration;
    
%%%%%% Load additive calibration
temp = load('camera_correcting/add_corr.mat');
add_corr = temp.add_corr;
 

rows = 32;
cols = 32;
for i_id = 1:length(ids)
    %%%%%% Load data
    datapath = ['/Data/' num2str(ids(i_id))];
    
    timing_response = double(h5read(file, datapath));
    wav = num2str(h5readatt(file, datapath, 'Wavelength [nm]'));
    source_no = h5readatt(file, datapath, 'Channel');
    rep_no = h5readatt(file, datapath, 'Run');
    target_pos = h5readatt(file, datapath, 'Targetted Head position {x, y, z, dist}');
    file_pos = ['target_' num2str(target_pos(1)) '_' ...
        num2str(target_pos(2)) '_' ...
        num2str(target_pos(3)) '_' ...
        num2str(target_pos(4)) ];
    if exist([fldr '/timing_data_corrected/' file_pos])==0
        mkdir([fldr '/timing_data_corrected/' file_pos])
    end
    %sig_interp = 100;
    %timing_response_1 = timing_response;
%     timing_data_tmp = zeros(size(timing_response)); % added by jingjing 2021.05.11
    for col_i = 1:cols
         
        for row_i = 1:rows

            if isempty(timing_response(:,row_i,col_i)) || (col_i == 32) &&...
                    (row_i == 2) || isnan(add_corr(row_i,col_i)) ||...
                   ((row_i==17&&col_i==24)) % changed by jingjing 2021.10.11
                timing_response(:, row_i, col_i) = 0; % added by jingjing 2021.05.11
            else
                timing_response(:, row_i, col_i) = correct_histogram(timing_response(:, row_i, col_i), add_corr(row_i, col_i), col_cal(col_i));
             
            end
            
        end
%         col_time = toc;
    end
     foo = sprintf([fldr '/' 'timing_data_corrected/' file_pos '/timing_response_%s_%s_%s_%s.mat'], ...
        Prefix, wav, num2str(source_no), num2str(rep_no));
%     disp(foo)
    save(foo, 'timing_response');
    %h5create(file, ['/Data_corrected/' num2str(ids(i_id))], [256 32 32], 'Datatype', 'uint32', 'ChunkSize', [256, 32, 32], 'Deflate', 2)
    %h5write(file, ['/Data_corrected/' num2str(ids(i_id))], timing_response)
end
toc
 
