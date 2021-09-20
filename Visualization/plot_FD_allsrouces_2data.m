%>@brief Brief description of the function
%>
%>visulize two Fourier domain data for all sources and detectors
%>
%>@param dat_1, dat_1 forward result
%>          .amplitude
%>          .phase
%>@param paras parameters  
%> 
%> author: jingjing jiang jing.jing.jiang@outlook.com
function plot_FD_allsrouces_2data(dat_1, dat_2, paras)
detnum = size(paras.det.coord,1);
srcnum = size(paras.src.coord,1);
figure;
for isrc = 1:srcnum
    id_dets = [1: detnum]+ (isrc-1)*detnum;
    subplot(2,srcnum,isrc)
    scatter3(paras.det.coord(:,1),paras.det.coord(:,2),...
        log(dat_1.amplitude(id_dets)) - mean(log(dat_1.amplitude(id_dets))), 'bo')
   hold on
   scatter3(paras.det.coord(:,1),paras.det.coord(:,2),...
        log(dat_2.amplitude(id_dets)) - mean(log(dat_2.amplitude(id_dets))), 'ro')
   
    title(['src' num2str(isrc) ' intensity'])
    subplot(2,srcnum,isrc + srcnum)
    scatter3(paras.det.coord(:,1),paras.det.coord(:,2),...
        dat_1.phase(id_dets) - mean(dat_1.phase(id_dets)),'o')
    hold on
    scatter3(paras.det.coord(:,1),paras.det.coord(:,2),...
        dat_2.phase(id_dets) - mean(dat_2.phase(id_dets)),'o')
    title('phase')
end