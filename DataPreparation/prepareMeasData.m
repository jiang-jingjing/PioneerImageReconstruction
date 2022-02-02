%>@brief Brief description of the function
%>
%> select time window and detectors
%>@param pathData path for reference measurements
%>@param paras parameters   
%>@param  bins array: selected time gates
%>@param cfg  config variable
%>@param isSavePosition save the selected positions  
%>@
%> @retval pos  positions of detectors updated
%> @retval data measured data
%> @retval paras parameters updated
%> author: jingjing jiang jing.jing.jiang@outlook.com
function [pos data paras] = prepareMeasData(pathData, paras, bins,...
    cfg, isSavePosition, varargin)

nPhoton_thresh = 0;

if ~isempty(varargin)
    if length(varargin) >= 1
        % log file name
        if isnumeric(varargin{1}) 
             nPhoton_thresh = varargin{1};
        else
            error('Bad 7th argument value. number expected for SmoothDelta.')
        end
    end
end
 % positions are saved in the calibration folder
fldr_cal = paras.calibration.dataPath;
fldr_pos = [fldr_cal{1}  '/DetSrcCoordinates/'];
fln_in = 'posSrcDet_Flat.mat';
load([fldr_pos '/' fln_in], 'pos')

% process the data from pathData
fldr = pathData{1};
flnm = {pathData{2}, pathData{3}};

len_bin = length(bins);
srcID =paras.src.num;
irep = paras.repetitionID;
detnum = size(pos.camera,1);
for  ii = 1:length(srcID)
   id_src =  srcID(ii);
    if irep>0 
        load([fldr flnm{1} num2str(id_src) '_' num2str(num_rep) flnm{2}]);

    else
         load([fldr flnm{1} num2str(id_src) flnm{2}]); % combined
    end

    xFOV = pos.camera(:,2); % flip x, y of the pixels
    yFOV = pos.camera(:,1); % changed by jingjing 2020.10.26
    hist_meas_temp = zeros(detnum,len_bin);
    count = 1;
    for jj = 1:length(xFOV)
        timing_tmp =  squeeze(timing_response(:, xFOV(jj), yFOV(jj)));
        if ~isempty(timing_tmp) & sum(timing_tmp)> nPhoton_thresh

            timing_tmp = timing_tmp - median(timing_tmp);% add back 2021.06.18
            timing_tmp_sft = timing_tmp(bins); % 2021.06.18
            hist_meas_temp(jj,:) =  timing_tmp_sft; 
            idValid(id_src).detid(count) =  jj ;
            count = count+1;
             hist_all(id_src).data = ...
                hist_meas_temp;
        end
    end
end

id_selected = idValid(1).detid;
srcnum = length(srcID);
if size(srcID,1)==1
    src_tmp = srcID;
elseif size(srcID,1)>1
    src_tmp = srcID';
end
for ii = src_tmp(2:end)
    id_selected = intersect(id_selected,idValid(ii).detid);
end

detnum = length(id_selected);
data.tpsf = zeros(detnum*srcnum, len_bin);
 
for isrc =src_tmp
    data.tpsf(1+(isrc-1)*detnum : isrc*detnum,:)= ...
     hist_all(isrc).data(id_selected,:);
%         idValid(id_src).detid= id_selected;
        
end       

pos.det.coordinates = pos.det.coordinates(id_selected, :);  
pos.camera = pos.camera(id_selected, :); 
 

if isSavePosition
% fln_out = 'detSrcPosFEM_selected.mat';
% save([fldr_pos '/' fln_out],  'PosSrcFEM', 'PosDetFEM',...
%     'pixel_size','IndexInFOVPiccolo', ...
%     'PosDetCam','id_selected');
% disp(['Selected FEM model positions are saved to ' [fldr_pos '/' fln_out]])
    fln_out = 'posSrcDet_Flat_selected.mat';
    save([fldr_pos '/' fln_out], 'pos')
    disp(['Selected FEM model positions are saved to ' [fldr_pos '/' fln_out]])

else
    disp(['do not save FEM model positions']);
end


data.time = cfg.tstart+cfg.tstep/2:cfg.tstep:cfg.tend;

end