%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo_nirfast_2d.m
% demo for nirfast simulation
% created on 2022.03.02 Jingjing Jiang
% modified on 2022.03.16 Jingjing Jiang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add paths for required packages
% define nirfast paths
% Nirfaster: GPU-fascilitated model
% Nirfast8: old CPU version
pathNirfaster = '/media/jiang/WD10T/Data//SoftwarePackages/NIRFASTer';
pathNirfast8 = '/media/jiang/WD10T/Data/SoftwarePackages/nirfast8';
addpath(genpath(pathNirfast8))
pathPioneerIR = '/media/jiang/WD10T/Data/Projects/PioneerImageReconstruction';
addpath(genpath(pathPioneerIR))
%% load mesh  / create mesh
mesh = load_mesh('./mesh/Rectangle-stnd-mesh');
mesh_homo = mesh;
%% remove short source-detector distances
nsrc = length(mesh_homo.source.num);
ndet = length(mesh_homo.meas.num);
psrc = mesh_homo.source.coord;
pdet = mesh_homo.meas.coord;
link_all_selected = ones(ndet*nsrc,1);
thresh_d = 10;
for isrc = 1:nsrc
    dist_mesh(isrc).value = ...
        sqrt(sum((repmat(psrc(isrc,:),ndet,1) - pdet).^2,2));
    idx_delete = find(dist_mesh(isrc).value<thresh_d);
    link_all_selected((isrc-1)*ndet+idx_delete) = 0;
end
id_sel = find(link_all_selected);
mesh_homo.link = mesh.link(id_sel,:);

%% visualize mesh

h_mesh = figure()
% plot boundary nodes
mesh_1 = mesh_homo; 
plot(mesh_1.nodes(logical(mesh_1.bndvtx),1),mesh_1.nodes(logical(mesh_1.bndvtx),2),'c.');
hold on
% plot sources
plot(mesh_1.source.coord(:,1),mesh_1.source.coord(:,2),'rx','LineWidth',2,'MarkerSize',8);
plot(mesh_1.source.coord(1,1),mesh_1.source.coord(1,2),'mx','LineWidth',2,'MarkerSize',12);
% plot detectors
plot(mesh_1.meas.coord(:,1), mesh_1.meas.coord(:,2),'bo','LineWidth',2,'MarkerSize',8);
plot(mesh_1.meas.coord(1,1), mesh_1.meas.coord(1,2),'co','LineWidth',2,'MarkerSize',10);

%% add optical properties
wavList = {'690', '830'};
 
% muas_bulk = [0.008 0.01];% Target values
% mus_r_bulk = [0.3 0.3]; % Target values
muas_bulk = [0.008 0.01]*2;% Target values
mus_r_bulk = [0.3 0.3]; % Target values
%  mus_r_bulk - 0.54;
muas_si = muas_bulk .*5;
mus_r_si = mus_r_bulk;
% mus_r_bulk = 0.54;
iwav = 1;
val.mua=muas_bulk(iwav); %692
val.mus=mus_r_bulk(iwav);
 
val.ri=1.43;
mesh_homo = set_mesh(mesh_homo,0,val);

%% add blob
clear blob
blob.x= 6;
blob.y= 5;
blob.r= 3;
blob.mua=muas_si(iwav);
blob.mus= mus_r_si(iwav);
blob.region=2;
mesh_anom = add_blob(mesh_homo,blob);

blob.x= -6;
blob.y= 5;
blob.r= 3;
blob.mua=muas_si(iwav);
blob.mus= mus_r_si(iwav);
blob.region=2;
mesh_anom = add_blob(mesh_anom,blob);

id_inc = find(mesh_anom.mua>0.008);
length(mesh_anom.nodes(id_inc,:))
plotimage(mesh_anom, mesh_anom.mua);

%% reduce mesh node size
mesh_homo.nodes = mesh_homo.nodes(:,1:2);
mesh_anom.nodes = mesh_anom.nodes(:,1:2);
%% calculate moments
addpath(genpath(pathNirfaster))
% test up to 2nd moment order (0th, 1st, and 2nd)
max_order = 2;

% increase the absolute tolerance for moments as the values will be low.
% The relative tolerance stays the same.
OPTIONS = solver_options;
OPTIONS.tolerance = 1e-12;

% GPU
tic;
[data_Moments] = femdata_stnd_TR_moments(mesh_homo, ...
    max_order,'field',[],solver_name_GPU,OPTIONS);
[data_Moments_anom] = femdata_stnd_TR_moments(mesh_anom, ...
    max_order,'field',[],solver_name_GPU,OPTIONS);
figure,
semilogy(data_Moments.moments)
hold on
semilogy(data_Moments_anom.moments)

%% reconstruction: moments
% reconstruction
tic
REG_THIKONOV = 1;
MESH_COARSE = [50 30]*0.5;
MAX_ITERATIONS = 10;
[meshRec_moments, pj_error_moments] = reconstruct_stnd_moments(mesh_homo, max_order, ...
     data_Moments_anom, 'mua',[], MESH_COARSE, MAX_ITERATIONS, REG_THIKONOV);
toc
plotimage(meshRec_moments, meshRec_moments.mua)

%% simulate frequency domain data 
% The unit of frequency in Nirfaster is Hz
% The unit of frequency in Nirfast8 is MHz
addpath(genpath(pathNirfaster))
frequencies = 100e6; % Hz
tic; [data_FD_NIRFAST_anom] = femdata_stnd_FD(mesh_anom,frequencies); 
figure,
subplot(121)
plot(data_FD_NIRFAST_anom.amplitude)
title('amplitude')
subplot(122)
plot(data_FD_NIRFAST_anom.phase)
title('phase')
%% reconstruction: frequency domain 
% data_FD_NIRFAST_anom = add_noise(data_FD_NIRFAST_anom, 2,2);
REG_THIKONOV =1;
tic
[meshRec_1F, pj_error_1F] = reconstruct_stnd_FD(mesh_homo,frequencies(1), ...
                             data_FD_NIRFAST_anom, 'mua',...
                            [],MESH_COARSE, MAX_ITERATIONS, REG_THIKONOV) 
toc
meshRec_1F
plotimage(mesh_anom, meshRec_1F.mua)

%% calculate temporal data
addpath(genpath(pathNirfaster))

% mesh_homo.link = mesh.link(find(link_all_selected),:);
% mesh_anom.link = mesh_homo.link;
cfg.tstart = 0;
cfg.tstep = 0.0488*1e-9 /2 ;  
cfg.num_gates = 20 *2;
cfg.tend = cfg.tstep* cfg.num_gates ;
cfg.timeGates = cfg.tstart+cfg.tstep:cfg.tstep:cfg.tend;
tic; [data_TR_anom] = femdata_stnd_TR(mesh_anom,cfg.tend,cfg.tstep,'field', 'BiCGStab_GPU'); toc;   
tic; [data_TR_homo] = femdata_stnd_TR(mesh_homo,cfg.tend,cfg.tstep,'field', 'BiCGStab_GPU'); toc;  
% tic; [data_TR_anom] = femdata_stnd_TR(mesh_anom,cfg.tend,cfg.tstep,'field', 'BiCGStab_CPU'); toc;   
% tic; [data_TR_homo] = femdata_stnd_TR(mesh_homo,cfg.tend,cfg.tstep,'field', 'BiCGStab_CPU'); toc; 
%%  plot tpsfs 
 itpsf = 6
h_psf=figure()
plot(cfg.timeGates, data_TR_homo.tpsf(itpsf,:))
hold on
plot(cfg.timeGates, data_TR_anom.tpsf(itpsf,:))

%% image reconstruction: Time domain
t_cfg.range = cfg.tend;
t_cfg.step = cfg.tstep;
t_cfg.num_gates = cfg.num_gates;
tic
[meshRec_TD, pj_error_TG] = reconstruct_stnd_TD(mesh_homo,t_cfg, ...
        data_TR_anom, 'mua',[], MESH_COARSE, MAX_ITERATIONS, REG_THIKONOV)
     toc
    plotimage(meshRec_TD, meshRec_TD.mua)

