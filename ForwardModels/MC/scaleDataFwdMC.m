%>@brief Brief description of the function
%>
%> generate tpsf for new absorption properties
%> based on monte carlo simulation
%>
%>@param resultMC 
%>@param paras parameters  
%>@param cfg  struct config parameter for time gates 
%>          .tstart
%>          .tend
%>          .tstep 
%>          .prop optical properties
%> 
%> @retval dataFwd forward result
function dataFwd = scaleDataFwdMC(resultMC, paras, cfg)
c0=299792458000 ; % [mm/s]
n = cfg.prop(2:end,4);
srcnum = size(paras.src.coord,1);
detnum = size(paras.det.coord,1);
dataFwd.time = cfg.tstart+cfg.tstep/2:cfg.tstep:cfg.tend;
len_bin = length(dataFwd.time);

for isrc = 1:srcnum


    detp =  resultMC(isrc).detp;
  
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
