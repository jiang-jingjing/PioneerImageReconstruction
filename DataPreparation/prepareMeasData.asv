%>@brief Brief description of the function
%>
%> get positions of sources and detectors from measruement
%>
%>@param volIn tissue volume
%>@param paras parameters           
%>@param cfg  bins array: selected time gates
%>
%> @retval pos  positions of detectors updated
%> @retval data measured data
%> @retval paras parameters updated

function [pos data paras] = prepareMeasData(pathData, paras, bins)

% positions are saved in the calibration folder
fldr_cal = paras.calibration.dataPath;
fldr_pos = [fldr_cal  '/DetSrcCoordinates/'];
fln_in = 'posSrcDet_Flat.mat';
load([fldr_pos '/' fln_in], 'pos')

% process the data from pathData
fldr = pathData{1};
flnm = {pathData{2}, pathData{3}};

idValid 
bins

pos.det.coordinates = pos.det.coordinates(idValid, :);  
pos.det.coordinates = pos.camera(idValid, :); 
 
fln_out = 'posSrcDet_Flat_selected.mat';
save([fldr_pos '/' fln_out], 'pos')

end