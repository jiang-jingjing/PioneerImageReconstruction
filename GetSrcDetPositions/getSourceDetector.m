%>@brief Brief description of the function
%>
%> get positions of sources and detectors from measruement
%>
%>@param volIn tissue volume
%>@param paras parameters           
%>
%>
%> @retval pos positions of detectors
%> @retval paras parameters updated
%> author: jingjing jiang jing.jing.jiang@outlook.com
function [pos paras] = getSourceDetector(pathData, paras)
fldr = pathData{1};
flnm = {pathData{2}, pathData{3}};
fldr_pos = [fldr  '/DetSrcCoordinates/'];
pixel_size = paras.det.pixel;
ratioR = paras.det.ratioR;
fovC = paras.det.fovC;
srcList = paras.src.num;
modelCenter = paras.det.modelCenter;
    if strcmp(paras.probe.type, 'Flat') || ...
            strcmp(paras.probe.type, 'flat')
        pos = GetSrcDetPositions_flat(fldr,  flnm, fldr_pos,...
            pixel_size, srcList,...
            -1, ratioR, fovC, modelCenter)
    else
        error('code for get coordinates of sources and detectors for flexible probe is not supported')
        
    end
end



