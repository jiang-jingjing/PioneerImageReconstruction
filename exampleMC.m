%>@brief Brief description of the function
%>
%> An example of image reconstruction for time domain modality 
%> based on MC
%> author: jingjing jiang jing.jing.jiang@outlook.com

%% ADD PATHS (temporary)
addpath(genpath('./'))
path_mcxlab = '/media/jiang/WD10T/Data/SoftwarePackages/mcx2021/';
addpath(genpath(path_mcxlab))
% mc 
% fld_my = '/media/jiang/WD10T/Data/Projects/NIROTreconstruction/RECONSTRUCTION_MATLAB';
% addpath([fld_my '/src/']);
% addpath([fld_my '/config/']);

% nirfast
% addpath('/media/jiang/WD10T/Data/Projects/NIROTreconstruction/RECONSTRUCTION_NIRFAST2019/')

%% flags
isCombRep = 0
isGenerateMesh = 1

%% folders
FOLDER = ['/media/jiang/WD10T/Data/Projects/NIROTreconstruction/PioneerMeasurement2021/Data/']
flag_prj = '/lp_20210615_depth2incl/lq_dist_11_6/';
fldr = [FOLDER flag_prj];
% file = [fldr '/laser-lqp_longexp.hdf5'];
Prefix = 'laser';
waveList = [689 725  ];
% ISS measured result
muas_bulk = [0.004  0.0045   ]  ;
mus_r_bulk = [0.648  0.593  ] ;  % estimated 2021.06.22 [mm-1]

srcList =  [1:11];
depth_list = [10  15 20 25 30]; % depth [mm]
distance_list = 8;%  lateral distance [mm]
rep_list = 1:100;
group_list = {'target_48.5_-1.5_26_-41.76',... % depth > 71 mm
   'target_48.5_66_26_-41.76',... % depth = 5 mm
   'target_48.5_61_26_-41.76',... % depth = 10 mm
   'target_48.5_56_26_-41.76',...% depth = 15 mm
   'target_48.5_51_26_-41.76',...% depth = 20 mm
   'target_48.5_46_26_-41.76'}; % depth = 25 mm

%% STAGE 1 Preparation 
%% Step 1a: define imaging modality 
nirot.ImagingModality = 'TD';
nirot.ReconstructionDataType = 'FD'; % single frequency
%% Step 1b: define forward model
% add corresponding forward package
nirot.Model = 'MC';
%% Step 1c: define tissue  
% nirot.Volume.Name = 'colin27_v3';
% nirot.Volume.Path = '/media/jiang/WD10T/Data/SoftwarePackages/mcx2021/mcxlab/examples/colin27_v3.mat';
nirot.Volume.Name = 'Cylinder';
nirot.Volume.Path = [];
%% Step 1d: get positions of sources and detectors
% get the positions from the Martin tapping experiment
% (homogeneous tissue)
nirot.probe.type = 'Flat';
nirot.det.pixel = 1.09;
nirot.det.fovC = [15 15] % approximate center of FOV
nirot.det.ratioR = 0.87; % ratio of FOV region
% nirot.det.modelCenter = [45 45];
nirot.det.modelCenter = [30 30];

nirot.wavelengths = [689 725 ];
nirot.iwav = 2;
% measured data for calibration
i_group = 1;
wav = num2str(nirot.wavelengths(nirot.iwav))
group = group_list{i_group};
fldr_r = [fldr '/' 'timing_data_corrected/' group]
flnm_h_r={['/timing_response_' Prefix '_' ...
    wav '_'], '.mat'};
nirot.calibration.dataPath = [fldr_r flnm_h_r];
nirot.src.num = [1:11];
nirot.repetitionID = 0; % no repetitions
[pos2D nirot] = getSourceDetector(nirot.calibration.dataPath,...
    nirot);

%% Step 1e: prepare measured data and select detectors
% specify time gates
len_bin =  50
bin0 = 13;
bin_select = bin0:bin0+len_bin-1; 
isSavePosition = 1;

% specify time
cfg.tstart=0;
cfg.tstep=0.0488e-9;
cfg.tend=cfg.tstep*len_bin;%5e-09;

nPhoton_thresh = 1e2;

[pos2D dataRef nirot] = prepareMeasData(...
    nirot.calibration.dataPath,...
    nirot,...
    bin_select,...
    cfg, ...
    isSavePosition,...
    nPhoton_thresh);
nirot.calibration.data = dataRef;
 
%% STAGE 2: Forward simulation
%% Step 2a: create tissue volume / mesh
filename_vol = ['./example_2inc/'  nirot.Volume.Name '.mat']
nirot.unitmm = 1;
nirot.vol = createCylinderMC(60,30,nirot.unitmm , filename_vol);

%load Colin27 brain atlas
% load nirot.Volume.Path
% nirot.vol = colin27

%% Step 2b: add optical properties
filename_vol_tiss = ['./example_2inc/'  nirot.Volume.Name ...
    '_opticalProperties.txt']
g = 0.35;
% g=0.01
mus  =  mus_r_bulk ./ (1-g);
n = 1.37;
% nirot.prop=[         0         0    1.0000    1.0000 ;% background/air
%     muas_bulk(nirot.iwav)    mus(nirot.iwav)   g    n % liquid
%     ];
nirot.prop=[         0         0    1.0000    1.0000 ;% background/air
    1e-6    mus(nirot.iwav)   g    n % liquid
    ];
fileID = fopen(filename_vol_tiss, 'w');
fprintf(fileID,'%1.4f %1.4f %1.2f %1.2f\n',nirot.prop');
fclose(fileID);

%% Step 2c: add sources and detectors to the volume surface
%load source detector from Step 1d
nirot = addSrcDetMC_Flat(pos2D, nirot);
% plot volume and src/det
h_vol = plotVolMC(nirot.vol, nirot);


% convert time domain data to Fourier domain
frequencies = 100 * 1e6; % Hz
tic
dataRef = td2fd(dataRef, frequencies);
toc
% visualie 2D distribution of log amplitude and phase
plot_FD_allsrouces(dataRef, nirot)
 
%% Step 2d: calculation of forward results
% MCX simulation
vol_init = nirot.vol;
cfg.nphoton=1e8;
cfg.maxdetphoton = 1e8;

% forward 
tic
[dataBase, nirot, cfg, resultMC] = forwardTimeMC(vol_init, ...
    nirot, cfg);
toc

h_fwd = figure; 
semilogy(dataBase.tpsf')

% frequencies = 100 * 1e6; % Hz
% tic
% dataBase= td2fd(dataBase, frequencies);
% toc
% plot_FD_allsrouces_2data(dataRef, dataBase, nirot)

%% generate tpsf for new absorption properties
% nirot.prop=[         0         0    1.0000    1.0000 ;% background/air
%     muas_bulk(nirot.iwav)    mus(nirot.iwav)   g    n % liquid
%     ];
cfg_2 = cfg;
cfg_2.prop(2,1) = muas_bulk(nirot.iwav);
tic
[dataFwd_2] = scaleDataFwdMC(resultMC, nirot, cfg_2);
toc
figure(h_fwd)
hold on
semilogy(dataFwd_2.tpsf')
 
%% convert time domain data to Fourier domain
frequencies = 100 * 1e6; % Hz
tic
dataFwd_2 = td2fd(dataFwd_2, frequencies);
toc
% visualie 2D distribution of log amplitude and phase
% plot_FD_allsrouces(dataFwd_2, nirot)
 
plot_FD_allsrouces_2data(dataRef, dataFwd_2, nirot)

%% convert time domain data to moments
n_moments = 3;
dataFwd = td2moments(dataFwd_2, n_moments);
figure, semilogy(dataFwd.moments )
% %% plotting TPSF at a voxel  
% figure
% % detpos=[30 14 9]; % choose a point inside the domain
% isrc = 4
% detid = 191;
% detpos =floor(nirot.det.coord(detid,:));
% detnum = size(nirot.det.coord,1);
% srcpos = nirot.src.coord(isrc,:);
% srcpos
% detpos
% Reff = 0.45 ; % to be calculated
% c0=299792458000 ; % [mm/s]
% n = nirot.prop(2:end,4);
% 
% f1 = resultMC(isrc).flux.data;
% twin=cfg.tstart+cfg.tstep/2:cfg.tstep:cfg.tend;
% tpsf_mc = squeeze(f1(detpos(1),detpos(2),detpos(3),:));
%  tpsf_de = tddiffusion(nirot.prop(2,1), ...
%     nirot.prop(2,2)*(1-nirot.prop(2,3)), ...
%     v, Reff, srcpos, detpos+1, twin);
% 
% % plot(twin*1e9,tpsf_de ./max(tpsf_de),'r');
% % hold on
% % plot(twin*1e9,tpsf_mc ./max(tpsf_mc) ,'o');
% % plot(twin*1e9,tpsf_mc_my ./max(tpsf_mc_my) ,'og');
% % 
% plot(twin*1e9,tpsf_de,'r');
% hold on
% plot(twin*1e9,tpsf_mc,'o');
%  
% 
% 
% 
% 
% 
% 
% %%
% c0=299792458000 ; % [mm/s]
% v = c0/n;
% srcpos = nirot.src.coord(isrc,:);
% % detpos=[30 14 9]; % choose a point inside the domain
% Reff = 0.45 ; % to be calculated
% f1 = dataFwd.flux.data;
% twin=cfg.tstart+cfg.tstep/2:cfg.tstep:cfg.tend;
% f_mc = zeros(detnum, len_bin);
% f_de = f_mc;
%  tic 
% for detid = 1:detnum
%     detpos =floor(nirot.det.coord(detid,:));
% 
% % srcpos
% % detpos
% tpsf_mc = squeeze(f1(detpos(1),detpos(2),detpos(3),:));
% f_mc(detid,:) = tpsf_mc;
% tpsf_de = tddiffusion(nirot.prop(2,1), ...
%     nirot.prop(2,2)*(1-nirot.prop(2,3)), ...
%     v, Reff, srcpos, detpos+1, twin);
% f_de(detid,:) = tpsf_de;
% end
% toc
%  
% %%
% 
% cw_f_mc = sum(f_mc,2);
% cw_f_de = sum(f_de,2);
% figure,
% plot3(nirot.det.coord(:,1), ...
%     nirot.det.coord(:,2), ...
%     log(cw_f_mc./mean(cw_f_mc)),...
%     'ob')
% hold on
% plot3(nirot.det.coord(:,1), ...
%     nirot.det.coord(:,2), ...
%     log(cw_f_de./mean(cw_f_de)),...
%     'or')
% % plot3(nirot.det.coord(:,1), ...
% %     nirot.det.coord(:,2), ...
% %     log(cw./mean(cw)),...
% %     'og')
% 
% plot3(nirot.det.coord(:,1), ...
%     nirot.det.coord(:,2), ...
%     log(cw_w1 ./mean(cw_w1)),...
%     'ok')
% xlabel('x')
% xlabel('y')
%% STAGE 3: Image reconstruction
%% Step 3a: calibration of measured data

%% Step 3b: image reconstruction
nirot_2 = nirot;
nirot_2.prop = cfg_2.prop;
cfg_2.replaydet=0;  % replay all det and sum all
% cfg_2.replaydet=1;  % replay only the 2nd detector
%cfg_2.replaydet=3;  % replay only the 3rd detector
% cfg_2.replaydet=-1; % replay all det and save all
% test jacobian
[jac, varargout] = jacobianTimeMC(nirot_2.vol, nirot_2,...
    resultMC, cfg_2 )

% reconstructionFD


 