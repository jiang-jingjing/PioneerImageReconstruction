%> demo MC simulation using mcxlab
%>
%> created by jingjing jiang jing.jing.jiang@outlook.com
%>
%% STAGE 1 Preparation 
%% Step 1a: define imaging modality 
nirot.ImagingModality = 'TD';
nirot.ReconstructionDataType = 'FD'; % single frequency
%% Step 1b: define forward model
% add corresponding forward package
nirot.Model = 'MC';
nirot.vol=uint8(ones(60,60,60));
nirot.src.coord=[30 30 60];
nirot.prop=[0 0 1 1;0.005 1 0 1.37];
nirot.det.coord=[15 30 60 ];
nirot.det.pixel = 4;
nirot.unitmm = 1;
%% calculation of forward results
% MCX simulation
vol_init = nirot.vol;
cfg.nphoton=1e8;
cfg.maxdetphoton = 1e8;
cfg.tstart=0;
cfg.tend=3e-9;
cfg.tstep=5e-11;
% forward 
tic
[dataBase, nirot, cfg, resultMC] = forwardTimeMC(vol_init, ...
    nirot, cfg);
toc

h_fwd = figure; 
semilogy(dataBase.time, dataBase.tpsf')

 
%% generate tpsf for new absorption properties
cfg_2 = cfg;
cfg_2.prop(2,1) = 0.01;
tic
[dataFwd_2] = scaleDataFwdMC(resultMC, nirot, cfg_2);
toc
figure(h_fwd)
hold on
semilogy(dataBase.time, dataFwd_2.tpsf')