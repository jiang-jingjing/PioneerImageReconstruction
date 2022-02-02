function [CS_ori CS_rec]=plot_truthVSrecon(meshSeg_truth, meshSeg_recon, ...
    coords)
addpath('/media/jiang/WD10T/Data/Projects/Projects/GAN/GAN_NIROT/NIRFAST_DATA/cs_plotmesh 091/cs_plotmesh 091')
x = coords(1);
y = coords(2);
z = round(coords(3));

CS_ori = cs_plotmesh_new(meshSeg_truth, [x y z], 1)
 
CS_rec = cs_plotmesh_new(meshSeg_recon,[x y z], 1)
% h_rec = figure;

figure
subplot(221)
plot(squeeze(CS_ori.XY.op(round(x)+3,:,1)))
hold on
plot(CS_rec.XY.op(round(x)+3,:,1))
legend('original','reconstructed')
title('dim 1') 
subplot(222)
plot(CS_ori.XZ.op(:,round(y)+3,1))
hold on
plot(CS_rec.XZ.op(:,round(y)+3,1))
title('dim 2') 
subplot(223)
cla
plot(CS_ori.XY.op(:,round(y)+3,1))
hold on
plot(CS_rec.XY.op(:,round(y)+3,1))
title('dim 1') 