%>@brief Brief description of the function
%>
%>visulie Fourier domain data for all sources and detectors
%>
%>@param dataFwd forward result
%>          .amplitude
%>          .phase
%>@param paras parameters  
%> 
function plot_FD_allsrouces(dataFwd, paras)
detnum = size(paras.det.coord,1);
srcnum = size(paras.src.coord,1);
figure;
for isrc = 1:srcnum
    id_dets = [1: detnum]+ (isrc-1)*detnum;
    subplot(2,srcnum,isrc)
    scatter3(paras.det.coord(:,1),paras.det.coord(:,2),...
        log(dataFwd.amplitude(id_dets)), 'o')
    title(['src' num2str(isrc) ' intensity'])
    subplot(2,srcnum,isrc + srcnum)
        scatter3(paras.det.coord(:,1),paras.det.coord(:,2),...
            dataFwd.phase(id_dets),'o')
   title('phase')
end