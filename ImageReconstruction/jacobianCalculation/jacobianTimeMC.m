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
%> author: jingjing jiang jing.jing.jiang@outlook.com 
function [jac, varargout] = jacobianTimeMC(volIn, paras, resultMC,cfg, varargin)

% % cfg.replaydet=0;  % replay all det and sum all
% % cfg.replaydet=1;  % replay only the 2nd detector
% %cfg.replaydet=3;  % replay only the 3rd detector
% cfg.replaydet=-1; % replay all det and save all

newcfg=cfg;
newcfg.outputtype='jacobian';
% numdet = size(paras.det.coord,1);
tic
for isrc = paras.src.num
newcfg.seed=resultMC(isrc).seeds.data;
 
newcfg.detphotons=resultMC(isrc).detp.data;
newcfg.srcpos = paras.src.coord(isrc,:);
% for idet = 1:numdet
%     cfg.replaydet=idet;
% cfg.replaydet = 2
[flux, detp, vol, seeds]=mcxlab(newcfg);
% 
% end
% J_tmp = flux.data;
% [v1 v2 v3 nt] = size(J_tmp); 
% nv = v1*v2*v3;
% J = reshape(J_tmp, nv, nt);
% jac =log10(abs(sum( reshape(J, v1,v2,v3,nt))));
jac=log10(abs(sum(flux.data,4)));
jac(~isfinite(jac)) = nan;


toc

figure, 
for iz = size(newcfg.vol,3) :-1:1
    imagesc(squeeze(jac(:,:,iz )))
    title(['ch ' num2str(isrc) ', z = ' num2str(iz)])
    colorbar
    pause(0.3)
end

figure
iy = round(size(newcfg.vol,1) /2);
imagesc(squeeze(jac(:,iy,:)))
end
end