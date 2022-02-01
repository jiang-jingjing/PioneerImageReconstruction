%>@brief Brief description of the function
%>
%> An example: prepare tpsf data for Pioneer
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
nirot.Volume.Name = 'Cylinder';
nirot.Volume.Path = [];
%% Step 1d: get positions of sources and detectors
% get the positions from the Martin tapping experiment
% (homogeneous tissue)
nirot.probe.type = 'Flat';
nirot.det.pixel = 1.09;%[mm]
nirot.det.fovC = [15 15] %[pixel] approximate center of FOV
nirot.det.ratioR = 0.87; % ratio of FOV region
% nirot.det.modelCenter = [45 45];
nirot.det.modelCenter = [30 30];

nirot.wavelengths = [689 725];
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
len_bin =  150
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
figure, plot(dataRef.time, dataRef.tpsf)
%% STAGE 2: Forward simulation
%% Step 2a: create tissue volume / mesh
filename_vol = ['./example_2inc/'  nirot.Volume.Name '.mat']
nirot.unitmm = 1;
nirot.vol = createCylinderMC(60,30,nirot.unitmm , filename_vol);


%% Step 2b: add optical properties
filename_vol_tiss = ['./example_2inc/'  nirot.Volume.Name ...
    '_opticalProperties.txt']
g = 0.35;
mus  =  mus_r_bulk ./ (1-g);
n = 1.37;
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