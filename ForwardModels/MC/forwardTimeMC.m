%>@brief Brief description of the function
%>
%> calculate forward result for time domain modality 
%> based on monte carlo simulation
%>
%>@param volIn tissue volume
%>@param paras parameters  
%>@param cfg  struct config parameter for time gates 
%>          .tstart
%>          .tend
%>          .tstep 
%> 
%> @retval dataFwd forward result
%> @retval paras parameters updated
%> @retval cfg MC config parameters  updated
%>
%> created by jingjing jiang jing.jing.jiang@outlook.com
%>
function [dataFwd, paras, cfg, varargout] = forwardTimeMC(volIn, paras,cfg )
 % shuffle
rng('shuffle')
cfg.seed= hex2dec('623F9A9E') * rand(1);
% volume
cfg.vol = volIn;
cfg.prop = paras.prop;

% sources
% % cfg.srctype='isotropic';
cfg.srcdir=[0 0 -1];
cfg.outputtype = 'flux';
% detectors
cfg.detpos = paras.det.coord;
cfg.detpos(:,4) =  paras.det.pixel ./ 2;
cfg.unitinmm = paras.unitmm;
% gpu specified
cfg.gpuid=1;
% cfg.gpuid='11'; % use two GPUs together
cfg.autopilot=1;
cfg.issrcfrom0=1;

cfg.isreflect=1; % enable reflection at exterior boundary
srcnum = size(paras.src.coord,1);
c0=299792458000 ; % [mm/s]
n = cfg.prop(2:end,4);
% v = c0./n;
detnum = size(paras.det.coord,1);
dataFwd.time = cfg.tstart+cfg.tstep/2:cfg.tstep:cfg.tend;
len_bin = length(dataFwd.time);
dataFwd.tpsf = zeros(srcnum*detnum, len_bin);
for isrc = 1:srcnum
    cfg.srcpos = paras.src.coord(isrc,:);

    [flux, detp, vol, seeds]=mcxlab(cfg);
    if  max(nargout,1)- 3>0

        resultMC(isrc).flux = flux;
        resultMC(isrc).detp = detp;
        resultMC(isrc).vol = vol;
        resultMC(isrc).seeds = seeds;
        clear flux detp vol seeds
    end
    % calculate tpsf
    times_seg = detp.ppath;
    detID = single(detp.detid);
    distancw_total=sum(detp.ppath,2) .* paras.unitmm;  
    bin = ceil(distancw_total  / (c0 ./n) / cfg.tstep); 
    bin(bin>len_bin) = len_bin;
    len = length(bin);
%     Track_atten=exp(-times_seg * muas_bulk(nirot.iwav)   * nirot.unitmm); % changed on 2018.03.09 
    Track_atten=exp(-times_seg * cfg.prop(2:end,1) * cfg.unitinmm);
    temp  = mexloop(len ,detID,bin,Track_atten);
    dataFwd.tpsf(1+(isrc-1)*detnum : isrc*detnum,:)= ...
        temp(1:detnum,1:len_bin); 
     
end
if max(nargout,1)-3>0
    varargout{1} = resultMC;
end
end