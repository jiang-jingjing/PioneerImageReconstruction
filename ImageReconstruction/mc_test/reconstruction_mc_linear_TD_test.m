%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image reconstruction based on Perturbasion model
% Forward model: MCXLAB time domain
% Inverse problem solver: linear, lsqr from Matlab
% Created by jingjing jiang, 
%            jing.jing.jiang@outlook.com
% Date: 2025.04.23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add paths
path_mcxlab = '/media/jiang/WD10T/Data/SoftwarePackages/mcx2023/';
addpath(genpath(path_mcxlab))
 

%% create an homogeneous case
%% Step 1 - run baseline simulation to get detected photons (and their seeds) at detectors
mua_homo = 0.005 * ones(1,2);
mua_heter = mua_homo;
mua_heter(2) = 0.1;

clear cfg
% cfg.seed = hex2dec('623F9A9E');
cfg.nphoton=1e8;
 
cfg.vol=uint8(ones(20,20,10));
cfg.vol(7:9, 11:13, 5:7)=2;
 

% time-domain simulation parameters
cfg.tstart = 0;
cfg.tend =0.3e-9;
cfg.tstep = 1e-11 ;
 

twin = cfg.tstart + cfg.tstep :cfg.tstep:cfg.tend;
numgate = length(twin);
 
cfg.prop=[0 0 1 1; 
    1e-6 0.4 0.5 1.37; 
    1e-6 0.4 0.5 1.37];
 
cfg.unitinmm = 1;

% add 10 x 10 detectors on the opposite side with the pixel size = 2 mm
% cfg.detpos=[30,20,0,1];        % to detect photons, one must first define detectors
det_spacing = 2;
 
x_start = 5; x_end = 15; % Includes center at 15
y_start = 5; y_end = 15;
z_pos =  10+1; % Opposite side of source (z=1)
XX = x_start:det_spacing:x_end;
YY = y_start:det_spacing:y_end;

% Generate detector coordinates
[x_det, y_det] = meshgrid(XX, YY);
cfg.detpos = [x_det(:), y_det(:), repmat(z_pos, numel(x_det), 1), ones(numel(x_det),1)];

% cfg.detpos = [10, 10,  10, 2 ];
cfg.maxdetphoton = 1e8;
cfg.issavedet =1;
 
cfg.autopilot = 1;
cfg.gpuid = 1;

cfg.isreflect = 0; % enable reflection at exterior boundaries
cfg.outputtype = 'flux'; 
 

%% run baseline simulation, now seed data are stored in output seeds.data
% srcList = [15, 15;
%     15 12 
%     15 18
%     12 15 
%     18 15] -5;
srcList = [10 10]
srcList = srcList  ;
srcList(:,3) =1;
num_src = size(srcList,1);
numdet = size(cfg.detpos,1);

%% Visualization Geometry
figure;
vol = cfg.vol;
slice(double(vol), 10, 10, 5); % Slices at x=10, y=10, z=5
shading interp
colormap(jet)
colorbar
axis tight
xlabel('X'), ylabel('Y'), zlabel('Z')
title('3D Volume Slices')
hold on
% Visualize detector placement (optional)
scatter3(cfg.detpos(:,1), cfg.detpos(:,2), cfg.detpos(:,3), 40, 'filled');
title('Detector Array (2 mm spacing)');
xlabel('X (voxels)'); ylabel('Y (voxels)'); zlabel('Z (voxels)');
axis equal; grid on;
% Add source position
hold on;
src_pos = srcList;
scatter3(src_pos(:,1), src_pos(:,2), src_pos(:,3), 150, 'filled', ...
    'MarkerFaceColor', [0 0.8 0], 'MarkerEdgeColor', 'k'); 
% axis equal tight
view(3)
%%
cw_1 = zeros(num_src*numdet,1);
cw_2 = cw_1;
for isrc = 1: num_src
% isrc = 1
cfg.srcpos = srcList(isrc,:) ;
cfg.srcdir = [0 0 1];
cfg.srctype='disk';
% cfg.srcparam1=[10, 0, 0, 0];
cfg.srcparam1=[5, 0, 0, 0]
[flux, detp, vol, seeds]=mcxlab(cfg);
% size(seeds.data)
seeds_all(isrc).seeds = seeds;
det_all(isrc).det = detp;
hist_1 = create_hist(detp, ...
    mua_homo, twin, cfg);
%% add anomaly  in the volume (heterogenous case)
 
hist_2 = create_hist(detp, ...
    mua_heter, twin, cfg);

%% 
numdet = size(cfg.detpos,1)
figure, 
for ii = 1:numdet
% ii = 4
    semilogy(hist_1(ii,:))
    hold on
    semilogy(hist_2(ii,:))
end

cw_1_tmp =  sum(hist_1,2);
cw_2_tmp =  sum(hist_2,2);

cw_1(1+numdet * (isrc-1):numdet * isrc) = cw_1_tmp;
cw_2(1+numdet * (isrc-1):numdet * isrc) = cw_2_tmp;

nx = length(XX);
ny = length(YY);
cw_2d_1 = reshape(cw_1_tmp,nx,ny);
cw_2d_2 = reshape(cw_2_tmp ,nx,ny);
figure,
subplot(131)
imagesc(cw_2d_1)
colorbar
subplot(132)
imagesc(cw_2d_2)
colorbar
subplot(133)
imagesc(cw_2d_1./cw_2d_2)
colorbar
end

 
 
%% Process measurements CW
% phi_homo = log10(cw_1);
% phi_hetero = log10(cw_2);
% data_diff = phi_hetero - phi_homo ;
%% Process measurements Time Gates
% remove zero gates
[idx_1 idx_2] = find(hist_1==0);
idx_start = max(idx_2) +1;
idx_tg_select = idx_start:numgate;
ng = length(idx_tg_select);
% [nd ng] = size(hist_1);

phi_homo = reshape(log10(hist_1(:,idx_tg_select))', numdet* ng , 1);
phi_hetero = reshape(log10(hist_2(:,idx_tg_select))', numdet* ng, 1);
data_diff = phi_hetero - phi_homo ;

 
% 
figure,
 
subplot(121)
plot(phi_homo )
hold on
plot(phi_hetero)

subplot(122)
plot(data_diff  )
%% Image reconstruction
 
%%  1. replay detected photon and output Jacobian [there is a bug with this replay mode
newcfg = cfg;
newcfg.prop(2:end,1) = mua_homo;
newcfg.outputtype = 'jacobian';
[vx, vy, vz] = size(newcfg.vol);
tic
jac_all = zeros(vx, vy, vz, numgate,numdet*num_src);
for isrc = 1:num_src
% newcfg.seed = seeds.data;
% newcfg.detphotons = detp.data;
% newcfg.replaydet = =-1 % 

newcfg.seed = seeds_all(isrc).seeds.data;
newcfg.detphotons = det_all(isrc).det.data;
jac = zeros(vx, vy, vz, numgate,numdet);
% for idet = 1: numdet
%     newcfg.replaydet= idet
    newcfg.replaydet= -1;
    [flux2, detp2, vol2, seeds2] = mcxlab(newcfg);
%     size(flux2.data)
    jac  = flux2.data;
    
% end
 jac_all(:,:,:,:,1+numdet * (isrc-1):numdet * isrc) = jac;
end
 
toc
%% 2. calculate jacobian matrix
% jac=log10(abs(jac1));
size(jac)

figure, 
% tg_plot = 10
for itg = 1:numgate
    imagesc(squeeze(jac(:,:,5,itg,1)))
    title(['time gate ' num2str(itg)])
    colorbar
    pause(0.3)
end
 
% jac(~isfinite(jac)) = nan;
%% reshape jacobian matrix
% J = reshape(jac, numel(cfg.vol), numdet)';  % [Ndet x Nvox]
% J = reshape(jac_all, numel(cfg.vol), numdet*num_src*numgate)' ;  % [Nmeas x Nvox]

% zero gates removed
J = reshape(jac_all(:,:,:,idx_tg_select,:), numel(cfg.vol), numdet*num_src*ng)' ;  % [Nmeas x Nvox]
  

size(J)

%% image reconstruction
 
tic
lambda = 1e0;  % Regularization parameter
%  JJt
H = J*J' + lambda*eye(size(J,1));
update = J'*(H\-data_diff);
update_3D = reshape(update, vx,vy,vz);
toc

 
%% jtj

% lambda = 1e-3;  % Regularization parameter
% H = J'*J + lambda*eye(size(J,2));
% update = H \ (J' * -data_diff);
% 
% update_3D = reshape(update, vx,vy,vz);

%% apply update and reconstruct
mu_a_recon_3D = update_3D+mua_homo(1);
figure
for iz = 1:size(cfg.vol,3)
    subplot(121)
    imagesc( mu_a_recon_3D(:,:,iz) )
    title(['z ' num2str(iz)])
    colorbar
    clim([0.005 nanmax(mu_a_recon_3D(:))])
    subplot(122)
    imagesc(cfg.vol(:,:, iz))

%     pause(0.3)
end
