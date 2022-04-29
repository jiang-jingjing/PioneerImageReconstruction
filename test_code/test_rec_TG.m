%% select detectors
addpath(genpath(pathNirfast8))

id_slct = mesh_homo.meas.num(1:10:end);
n_det = length(id_slct);
% mesh_2.meas.num = mesh_homo.meas.num(id_slct);
% mesh_2.meas.coord = mesh_homo.meas.coord(id_slct,:);

pos2D_2 = pos2D;
pos2D_2.det.coordinates = pos2D.det.coordinates(id_slct,:);

% link_all = [];
% for is = 1: 11
%     for id = 1: n_det
%         link_tmp = [is id 1];
%         link_all = [link_all;link_tmp];
%     end
% end
% mesh_2.link = link_all;
% mesh_2 = minband_opt(mesh_2);

% create mesh
mesh_2 = create_cylinder_D90(pos2D_2, 3, 1, 1); %D90 d50
% fn_mesh = 'test_mesh';
% save_mesh(mesh_2, fn_mesh);
% mesh_2 = load_mesh(fn_mesh);
%% - - Step 2d-2: datatype 2: time gate
cfg.tstart = 0;
cfg.tstep = 0.0488*1e-9;  
cfg.num_gates = len_bin;
cfg.tend = cfg.tstep* cfg.num_gates ;
cfg.timeGates = cfg.tstart+cfg.tstep:cfg.tstep:cfg.tend;

 
%% visualize mesh 3D
h_mesh = figure()
% plot boundary nodes
plot3(mesh_2.nodes(logical(mesh_2.bndvtx),1),...
    mesh_2.nodes(logical(mesh_2.bndvtx),2),...
    mesh_2.nodes(logical(mesh_2.bndvtx),3),'c.');
hold on
% plot sources
plot3(mesh_2.source.coord(:,1),mesh_2.source.coord(:,2),...
    mesh_2.source.coord(:,3),'rx','LineWidth',2,'MarkerSize',8);
% plot detectors
plot3(mesh_2.meas.coord(:,1), mesh_2.meas.coord(:,2),mesh_2.meas.coord(:,3),...
    'bo','LineWidth',2,'MarkerSize',8);
 
%% - Step 2c: add inclusions / or load segmentation 
%model with a spherical inclusion
mesh_2_heter = mesh_2;
xc = mean(mesh_2.nodes(:,1));
yc = mean(mesh_2.nodes(:,2));
zc = mean(mesh_2.nodes(:,3));
[xc yc zc]

R_inc = 5.1;
depth = 15;
% changed 2021.12.13
mask_temp1 = (mesh_2_heter.nodes(:,1)-xc).^2 +...
    (mesh_2_heter.nodes(:,3)-(max(mesh_2_heter.nodes(:,3))-depth)).^2 ...
    < R_inc*R_inc;
mask_temp2 = mesh_2_heter.nodes(:,2)<(yc+15) & ...
    mesh_2_heter.nodes(:,2)> (yc-15) ; 

mask_ano =  mask_temp1 .* mask_temp2;
 id_ano =  find(mask_ano) ;  %  30 x 10 
 length(id_ano)
 
mesh_2_heter.region(id_ano) = 2;

 
val_ano.mua= muas_inc(nirot.iwav);
val_ano.mus= mus_r_inc(nirot.iwav)% 
mesh_2_heter = set_mesh(mesh_2_heter,2,val_ano);

z0 = max(mesh_2_heter.nodes(:,3))- depth;
CS_ori = cs_plotmesh_new(mesh_2_heter, [xc yc z0], 1)
%% calculate temporal data

addpath(genpath(pathNirfaster))

tic; 
[dataSIM_2] = femdata_stnd_TR(mesh_2,...
    cfg.tend,cfg.tstep,'field', 'BiCGStab_GPU');
toc;   
%  plot TPSFs
h_tg=figure();
semilogy(cfg.timeGates, dataSIM_2.tpsf')
xlabel('time [s]')
title('TPSFs')

cfg.timeGates = cfg.tstart+cfg.tstep:cfg.tstep:cfg.tend;
[dataSIM_2_heter] = femdata_stnd_TR(mesh_2_heter,...
    cfg.tend,cfg.tstep,'field', 'BiCGStab_GPU');
%%  plot tpsfs 
itpsf = 6
h_psf=figure()
plot(cfg.timeGates, dataSIM_2.tpsf(itpsf,:))
hold on
plot(cfg.timeGates, dataSIM_2_heter.tpsf(itpsf,:))

%% image reconstruction: Time domain
t_cfg.range = cfg.tend;
t_cfg.step = cfg.tstep;
t_cfg.num_gates = cfg.num_gates;
MESH_COARSE = [90 90 50]/10;
MAX_ITERATIONS = 10;
REG_THIKONOV =10
tic
[meshRec_TD_SIM, pj_error_TG] = reconstruct_stnd_TD(mesh_2,t_cfg, ...
        dataSIM_2_heter, 'mua',[], MESH_COARSE, MAX_ITERATIONS, REG_THIKONOV)
     toc
 coords4plot = [xc yc  z0];
plot_truthVSrecon(mesh_2_heter, ...
            meshRec_TD_SIM,  coords4plot);    
