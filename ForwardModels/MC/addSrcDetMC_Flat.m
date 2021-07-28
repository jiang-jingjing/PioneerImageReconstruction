%>@brief Brief description of the function
%>
%> add sources and detectors to the volume surface 
%> for MC methods
%> for the flat probe
%> 
%>@param  pos2D sources and detectors 2D from measured data
%>@param  paras parameeters
%>
%> @retval paras update 3D coordinates of sources and detectors 
%>  in  parameters
function paras = addSrcDetMC_Flat(pos2D, paras)
    vol = paras.vol;
    [x y z] = size(vol);
    pos3D = pos2D;
    pos3D.det.coordinates(:,3) = z;
    pos3D.src.coordinates(:,3) = z;
    paras.det.coord = pos3D.det.coordinates;
    paras.src.coord = pos3D.src.coordinates;   
end

 