%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image reconstruction based on Perturbasion model
% Forward model: MCXLAB
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
cfg.nphoton=3e7;
cfg.vol=uint8(ones(20,20,10));
cfg.vol(10:12, 8:10, 3:5)=2;

% time-domain simulation parameters
cfg.tstart = 0;
cfg.tend =3e-9;
% cfg.tstep = 5e-11 ./5;
cfg.tstep = 3e-9;

twin = cfg.tstart + cfg.tstep :cfg.tstep:cfg.tend;
numgate = length(twin);

cfg.prop=[0 0 1 1; 
    1e-6 0.4 0.5 1.37; 
    1e-6 0.4 0.5 1.37];
 
cfg.unitinmm = 1;

det_spacing = 2;
x_start = 5; x_end = 15; % Includes center at 15
y_start = 5; y_end = 15;
z_pos =  1; % Opposite side of source (z=1)
XX = x_start:det_spacing:x_end;
YY = y_start:det_spacing:y_end;

% Generate detector coordinates
[x_det, y_det] = meshgrid(XX, YY);
cfg.detpos = [x_det(:), y_det(:), repmat(z_pos, numel(x_det), 1), ones(numel(x_det),1)];
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
    cfg.srcparam1=[1, 0, 0, 0]
    [flux, detp, vol, seeds]=mcxlab(cfg);
    % size(seeds.data)
    seeds_all(isrc).seeds = seeds;
    det_all(isrc).det = detp;
    hist_1 = create_hist(detp, ...
        mua_homo, twin, cfg);
    %% add anomaly  in the volume (heterogenous case)
    %  cfg_hetero = cfg;
    % cfg_hetero.prop(3,1) = 0.1;  % Add anomaly
    hist_2 = create_hist(detp, ...
        mua_heter, twin, cfg);
    
    %% prepare forward results
    % numdet = size(cfg.detpos,1)
    % figure, 
    % for ii = 1:numdet
    % % ii = 4
    %     semilogy(hist_1(ii,:))
    %     hold on
    %     semilogy(hist_2(ii,:))
    % end
    
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
phi_homo = log10(cw_1);
phi_hetero = log10(cw_2);
data_diff = phi_hetero - phi_homo ;
% %% Process measurements Time Gates
% phi_homo = log10(hist_1);
% phi_hetero = log10(hist_2);
% data_diff = phi_hetero - phi_homo ;
% 
figure,
 
subplot(121)
plot(phi_homo )
hold on
plot(phi_hetero)

subplot(122)
plot(data_diff' )

%%    replay detected photon and output Jacobian [there is a bug with this replay mode
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
tg_plot = 1
for idet = 1:numdet
    imagesc(squeeze(jac(:,:,10,tg_plot,idet)))
    title(['det ' num2str(idet)])
    colorbar
    pause(0.3)
end

%% reshape jacobian matrix

% J = reshape(jac, numel(cfg.vol), numdet)';  % [Ndet x Nvox]
J = reshape(jac_all, numel(cfg.vol), numdet*num_src')' ;  % [Nmeas x Nvox]
  
size(J)
 
%% Option 1: lsqr
 
reg_thikonov = 1e-3 ;
d=size(J);
L=eye(d(1,2),d(1,2));
tol = 1e-12;
 
maxit = 1e3

% 


[recon_update,fl,rr,it,rv,lsrv]=lsqr([...
    J;reg_thikonov*L],[-data_diff;zeros(d(1,2),1)],...
    tol,maxit,[],[],[]);


figure
semilogy(0:size(lsrv,1)-1,lsrv,'--o',0:size(rv,1)-1,rv,'-o')
legend("Least-squares residual","Relative residual")

recon_update_3D = reshape(recon_update, vx,vy,vz);
rec_3D = recon_update_3D+mua_homo(1);

figure
for iz = 1:size(cfg.vol,3)
    subplot(121)
    imagesc( rec_3D(:,:,iz) )
    title(['z ' num2str(iz)])
    colorbar
%     clim([0.003 0.4])
    subplot(122)
    imagesc(cfg.vol(:,:, iz))
    pause(0.3)
end
 
%% Option 2:  direct calculation
tic
lambda = 1e-3;  % Regularization parameter
%  JJt
H = J*J' + lambda*eye(size(J,1));
update = J'*(H\-data_diff);
update_3D = reshape(update, vx,vy,vz);
toc
%% jtj, when number of meas > number of voxels

% lambda = 1e-3;  % Regularization parameter
% H = J'*J + lambda*eye(size(J,2));
% update = H \ (J' * data_diff);

% update_3D = reshape(recon_update, vx,vy,vz);

%%   Apply update and reconstruct
mu_a_recon_3D = update_3D+mua_homo(1);
figure
for iz = 1:size(cfg.vol,3)
    subplot(121)
    imagesc( mu_a_recon_3D(:,:,iz) )
    title(['z ' num2str(iz)])
    colorbar
%     clim([0.005 0.015])
    subplot(122)
    imagesc(cfg.vol(:,:, iz))

    pause(0.3)
end

