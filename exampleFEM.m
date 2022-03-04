%>@brief Brief description of the function (not finished
%>
%> An example of image reconstruction for time domain modality 
%> based on FEM (Nirfast)
%> author: jingjing jiang jing.jing.jiang@outlook.com
%% Step 0
%% - add PATHS (temporary)
addpath(genpath('./'))
% define nirfast paths
% Nirfaster: GPU-fascilitated model
% Nirfast8: old CPU version
pathNirfaster = '/media/jiang/WD10T/Data//SoftwarePackages/NIRFASTer';
pathNirfast8 = '/media/jiang/WD10T/Data/SoftwarePackages/nirfast8';

%% - flags
isCombRep = 0; % 1: combine measured data

%% folders
FOLDER = ['/media/jiang/WD10T/Data/Projects/NIROTreconstruction/PioneerMeasurement2021/Data/']
Prefix = 'laser';
rep_list = 1:5;
group_list = {'target_47_71_19_0'}; 
depth_list = 15; %guess
i_group = 1;
group = group_list{i_group};
fldr_homo = [FOLDER   '/TM_SPih10r_ref/']
fldr_heter = [FOLDER   '/TM_SPih10_ref/']
%% optical properties 
% optical properties of the bulk
muas_bulk = 0.0055 * ones(1,4); %690 725 802 840
mus_r_bulk = [0.81  0.76 0.69 0.67] ;  % 
% optical properties of the inclusions
muas_inc = [0.0366 0.030 0.0187 0.124];
mus_r_inc = [0.87 0.9 0.831 0.8];
%% STAGE 1 Preparation 
%% - Step 1a: define imaging modality 
nirot.ImagingModality = 'TD';
% nirot.ReconstructionDataType = 'FD'; % single frequency
% nirot.ReconstructionDataType = 'TG'; % all time gates
nirot.ReconstructionDataType = 'Moments'; % temporal moments
%% - Step 1b: define forward model
% add corresponding forward package
nirot.Model = 'FEM';
%% - Step 1c: define tissue  
% nirot.Volume.Name = 'colin27_v3';
% nirot.Volume.Path = '/media/jiang/WD10T/Data/SoftwarePackages/mcx2021/mcxlab/examples/colin27_v3.mat';
nirot.Volume.Name = 'Cylinder';
nirot.Volume.Path = [];
%% - Step 1d: get positions of sources and detectors
% get the positions from the Martin tapping experiment
% (homogeneous tissue)
nirot.probe.type = 'Flat';
nirot.det.pixel = 1.09;
nirot.det.fovC = [17 15.5] % approximate center of FOV
nirot.det.ratioR = 0.87; % ratio of FOV region
nirot.det.modelCenter = [45 45];
% nirot.det.modelCenter = [30 30];

nirot.wavelengths = [690 725 802 839];
nirot.iwav = 2;
% measured data for calibration
wav = num2str(nirot.wavelengths(nirot.iwav))

fldr_r = [fldr_homo '/' 'timing_data_corrected/' group]
flnm_h_r={['/timing_response_' Prefix '_' ...
    wav '_'], '.mat'};
nirot.calibration.dataPath = [fldr_r flnm_h_r];
nirot.src.num = [1:11];
nirot.repetitionID = 0; % no repetitions
[pos2D nirot] = getSourceDetector(nirot.calibration.dataPath,...
    nirot);
% measured data for inclusion case
fldr_inc = [fldr_heter '/' 'timing_data_corrected/' group]
nirot.measurements.dataPath = [fldr_inc flnm_h_r];
%% - Step 1e: prepare measured data and select detectors
% specify time gates
len_bin =  50
bin0 = 13;
bin_select = bin0:bin0+len_bin-1; 
isSavePosition = 1;

cfg.tstart=0;
cfg.tstep=0.0488e-9;
cfg.tend=cfg.tstep*len_bin;%5e-09;

nPhoton_thresh = 1e2;

%% - - Step 1e-1: homogeneous case
[pos2D dataRef nirot] = prepareMeasData(...
    nirot.calibration.dataPath,...
    nirot,...
    bin_select,...
    cfg, ...
    isSavePosition,...
    nPhoton_thresh);
nirot.calibration.data = dataRef;
%% - - Step 1e-2: inclusion case
isSavePosition = 0;
[pos_tmp dataInc nirot] = prepareMeasData(...
    nirot.measurements.dataPath,...
    nirot,...
    bin_select,...
    cfg, ...
    isSavePosition,...
    nPhoton_thresh);
nirot.measurements.data = dataInc;
%% STAGE 2: Forward simulation
%% - Step 2a: create tissue volume / mesh with sources and detectors
% mesh for Nirfast
% generate / load mesh
isGenerateMesh = 1; 
if isGenerateMesh
    addpath(genpath(pathNirfast8))
    mesh = create_cylinder_D90(pos2D, 3, 1, 1); %D90 d50
    fn_mesh = 'exampleVolumeFEM/Cylinder_D90/Cylinder_1_mesh';
    nirot.Volume.Path = [fn_mesh  '.mat'];
    save(nirot.Volume.Path, 'mesh');

else % load mesh
    fn_mesh = 'exampleVolumeFEM/Cylinder_D90/Cylinder_1_mesh';
    nirot.Volume.Path = [fn_mesh  '.mat'];
    mesh=load(nirot.Volume.Path);
end
nirot.det = mesh.meas;
nirot.src = mesh.source;
%% - Step 2b: add optical properties to the mesh
g = 0.35;
mus  =  mus_r_bulk ./ (1-g);
n = 1.37;
val_p.mua=muas_bulk(nirot.iwav);  
val_p.mus=mus_r_bulk(nirot.iwav); 
val_p.ri=n; % pdms

n_region = unique(mesh.region);
mesh_homo = mesh;
for ii = n_region'
mesh_homo = set_mesh(mesh_homo,ii,val_p);
end

% backup the optical properties
nirot.prop=[  0  0 1.0000 1.0000; % background/air
    val_p.mua  mus(nirot.iwav)   g    n % liquid
    ];
filename_vol_tiss = ['./exampleVolumeFEM/'  nirot.Volume.Name ...
    '_opticalProperties.txt']
fileID = fopen(filename_vol_tiss, 'w');
fprintf(fileID,'%1.4f %1.4f %1.2f %1.2f\n',nirot.prop');
fclose(fileID);

%% - Step 2c: add inclusions / or load segmentation 
%model with a spherical inclusion
mesh_heter = mesh_homo;
xc = mean(mesh_homo.nodes(:,1));
yc = mean(mesh_homo.nodes(:,2));
zc = mean(mesh_homo.nodes(:,3));
[xc yc zc]

R_inc = 5.1;
depth = 15;
% changed 2021.12.13
mask_temp1 = (mesh_heter.nodes(:,1)-xc).^2 +...
    (mesh_heter.nodes(:,3)-(max(mesh_heter.nodes(:,3))-depth)).^2 ...
    < R_inc*R_inc;
mask_temp2 = mesh_heter.nodes(:,2)<(yc+15) & ...
    mesh_heter.nodes(:,2)> (yc-15) ; 

mask_ano =  mask_temp1 .* mask_temp2;
 id_ano =  find(mask_ano) ;  %  30 x 10 
 length(id_ano)
 
mesh_heter.region(id_ano) = 2;

 
val_ano.mua= muas_inc(nirot.iwav);
val_ano.mus= mus_r_inc(nirot.iwav)% 
mesh_heter = set_mesh(mesh_heter,2,val_ano);

z0 = max(mesh_heter.nodes(:,3))- depth;
CS_ori = cs_plotmesh_new(mesh_heter, [xc yc z0], 1)
%% - Step 2d: calculation of forward results
% add path for GPU-version Nirfaster
addpath(genpath(pathNirfaster))

%% - - Step 2d-1: datatype 1: Frequency domain
if strcmp(nirot.ReconstructionDataType, 'FD')
OPTIONS = solver_options;
OPTIONS.tolerance = 1e-12;
frequencies = 100*1e6;
tic; [dataSIM_fd_homo] = femdata_stnd_FD(mesh_homo,...
    frequencies(1), OPTIONS); 
toc; 
plot_data(dataSIM_fd_homo)
elseif strcmp(nirot.ReconstructionDataType, 'TG')
%% - - Step 2d-2: datatype 2: time gate
cfg.tstart = 0;
cfg.tstep = 0.0488*1e-9;  
cfg.num_gates = len_bin;
cfg.tend = cfg.tstep* cfg.num_gates ;
cfg.timeGates = cfg.tstart+cfg.tstep:cfg.tstep:cfg.tend;
tic; 
[dataSIM_tg_homo] = femdata_stnd_TR(mesh_homo,...
    cfg.tend,cfg.tstep,'field', 'BiCGStab_GPU');
toc;   
%  plot TPSFs
h_tg=figure();
semilogy(cfg.timeGates, dataSIM_tg_homo.tpsf')
xlabel('time [s]')
title('TPSFs')
elseif strcmp(nirot.ReconstructionDataType, 'Moments')
%% - - Step 2d-3: datatype 3: temporal moments
% simulation
max_order = 2;

% increase the absolute tolerance for moments as the values will be low.
% The relative tolerance stays the same.
OPTIONS = solver_options;
OPTIONS.tolerance = 1e-12;
tic
[dataSIM_m_homo] = femdata_stnd_TR_moments(mesh_homo, ...
    max_order,'field',[],solver_name_GPU,OPTIONS);
toc
% plot moments
h_m = figure;
semilogy(dataSIM_m_homo.moments)
title('moments')
end


%% STAGE 3: Image reconstruction
%% - Step 3a: datatype 1: Frequency domain
if strcmp(nirot.ReconstructionDataType, 'FD')
%% - - Step 3a-1: Convert TD to FD reference data
    tic
    dataEXP_fd_homo= td2fd(nirot.calibration.data, ...
        frequencies);
    toc
    % visualie 2D distribution of log amplitude and phase
    plot_FD_allsrouces(dataEXP_fd_homo, nirot)
%% - - Step 3a-2: calibration of FD data
    tic
    dataEXP_fd_heter= td2fd(nirot.measurements.data,...
        frequencies,...
        dataEXP_fd_homo,...
        dataSIM_fd_homo);
    toc    
%     plot_FD_allsrouces_2data(dataEXP_fd_homo, dataEXP_fd_heter, nirot)
    
%% - - Step 3a-3: simulation validation [FD]
    tic; [dataSIM_fd_heter] = femdata_stnd_FD(mesh_heter,...
    frequencies(1), OPTIONS); 
    toc; 
    mesh_coarse = [22 22 12.5];
    n_iter = 15;
    lambda =  10;
    % STND
    tic
    [meshRec_fd_SIM, pj_error] = reconstruct_stnd_FD(...
                mesh_homo,frequencies(1), ...
                dataSIM_fd_heter, 'mua',...
                [],mesh_coarse,n_iter, lambda);
    toc
    coords4plot = [xc yc  z0];
    plot_truthVSrecon(mesh_heter, ...
            meshRec_fd_SIM,  coords4plot); 
%% - - Step 3a-4: Image reconstruction [FD] measured
    mesh_coarse = [22 22 12.5];
    n_iter = 15;
    lambda =  10;
    % STND
    tic
    [meshRec_FD, pj_error] = reconstruct_stnd_FD(...
                mesh_homo,frequencies(1), ...
                dataEXP_fd_heter, 'mua',...
                [],mesh_coarse,n_iter, lambda);
    toc
    coords4plot = [xc yc  z0];
    plot_truthVSrecon(mesh_heter, ...
            meshRec_FD,  coords4plot); 
elseif strcmp(nirot.ReconstructionDataType, 'TG')
%% - Step 3b: datatype 2: time gate
%% - - Step 3b-1: calibration of TG data
%% - - Step 3b-2: simulation validation [TG]
%% - - Step 3b-3: Image reconstruction [TG]
elseif strcmp(nirot.ReconstructionDataType, 'Moments')
%% - Step 3c: datatype 3: moments
%% - - Step 3c-1: Convert TD to moments data
    n_moments = max_order+1;
    dataEXP_m_homo = td2moments(nirot.calibration.data, ...
        n_moments);
    figure, semilogy(dataEXP_m_homo.moments)
%% - - Step 3c-2: calibration of moments data
    dataEXP_m_heter = td2moments(nirot.measurements.data,...
        n_moments, dataEXP_m_homo, dataSIM_m_homo);
%% - - Step 3c-3: simulation validation [Moments]
    tic
    [dataSIM_m_heter] = femdata_stnd_TR_moments(mesh_heter, ...
        max_order,'field',[],solver_name_GPU,OPTIONS);
    toc
         
    figure
    semilogy(dataSIM_m_heter.moments)
    hold on
    semilogy(dataEXP_m_heter.moments)
    
    mesh_coarse = [22 22 12.5];
    n_iter = 15;
    lambda =  10;
    % STND
    [meshRec_m_SIM, pj_error] = reconstruct_stnd_moments(...
                mesh_homo,max_order, ...
                dataSIM_m_heter, 'mua',...
                [],mesh_coarse,n_iter, lambda);
    coords4plot = [xc yc  z0];
    plot_truthVSrecon(mesh_heter, ...
            meshRec_m_SIM,  coords4plot);       
%% - - Step 3c-4: Image reconstruction [Moments]
    % STND
    [meshRec_m_MEAS, pj_error] = reconstruct_stnd_moments(...
                mesh_homo,max_order, ...
                dataEXP_m_heter, 'mua',...
                [],mesh_coarse,n_iter, lambda);
    coords4plot = [xc yc  z0];
    plot_truthVSrecon(mesh_heter, ...
            meshRec_m_MEAS,  coords4plot);                
end