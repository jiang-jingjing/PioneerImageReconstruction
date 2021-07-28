%>@brief Brief description of the function
%>
%> calculate jacobian for time domain modality 
%> based on monte carlo simulation
%>
%>@param volIn tissue volume
%>@param paras parameters           
%>
%>
%> @retval jac jacobian
 
function [jac, varargout] = jacobianTimeMC(volIn, paras, resultMC,cfg, varargin)

% cfg.replaydet=0;  % replay all det and sum all
cfg.replaydet=1;  % replay only the 2nd detector
%cfg.replaydet=3;  % replay only the 3rd detector
% cfg.replaydet=-1; % replay all det and save all

newcfg=cfg;
newcfg.outputtype='jacobian';
for isrc = paras.src.num
newcfg.seed=resultMC(isrc).seeds.data;
 
newcfg.detphotons=resultMC(isrc).detp.data;
 
[flux, detp, vol, seeds]=mcxlab(newcfg);
jac=log10(abs(sum(flux.data,4)));
jac(~isfinite(jac)) = nan;

end


end