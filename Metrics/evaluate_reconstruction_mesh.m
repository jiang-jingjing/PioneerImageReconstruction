function [metric_eval mesh_segmented meshSeg_GT] = ...
    evaluate_reconstruction_mesh(...
    mesh_GT, mesh_rec,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return: metric_eval.rmse
%                    .psnr
%                    .dice
%                    .snr 
% SYNTAX:
% [metric_eval mesh_segmented meshSeg_GT] =
%   evaluate_reconstruction_mesh(mesh_GT, mesh_rec,[])
%	evaluate_reconstruction_mesh(mesh_GT, mesh_rec, thresh)
%   evaluate_reconstruction_mesh(mesh_GT, mesh_rec, thresh,'mua')
%   evaluate_reconstruction_mesh(mesh_GT, mesh_rec, thresh,'mua','resolution')
% 
% created 2021.05.10 by Jingjing Jiang jing.jing.jiang@outlook.com
% modified:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check in/out

narginchk(1,4);
nargoutchk(0,3);

% default

parameter = 'mua';
isResolution = 0; % does not check resolution,
thresh = nan;
mesh_segmented = [];
meshSeg_GT = [];
%1: calculate resolution (FWHM distance)
if ~isempty(varargin)
    if length(varargin) >= 1
        if ~isnumeric(varargin{1}) 
            % sanity check
            if varargin{1} < 0
                error('wrong parameter. Please see help for details on how to use this function.')
            else
                thresh = varargin{1};
            end
        else
%             error('Bad 3rd argument value. A string expected. Please see help for details on how to use this function.')
               disp('default thresholding')
               thresh = nan;
        end
    end
    if length(varargin) >= 2
        if ~ischar(varargin{2}) && ~isstring(varargin{2})
            % sanity check
            if varargin{1} < 0
                error('wrong parameter. Please see help for details on how to use this function.')
            else
                parameter = varargin{2};
            end
        else
            error('Bad 4th argument value. A string expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 3
        if ~ischar(varargin{3}) && ~isstring(varargin{3})
            % sanity check
            if varargin{3} < 0
                error('wrong input. Please see help for details on how to use this function.')
            else
                if strcmp(lower(varargin{3}), 'resolution')
                    isResolution = 1;
                end
            end
        else
            error('Bad 5th argument value. A string expected. Please see help for details on how to use this function.')
        end
    end    
    if length(varargin) > 4
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% calculate evaluation metrics
if strcmp(parameter,  'mua')
   disp( 'evaluate mua only')  
   MSE = mean((mesh_GT.mua - mesh_rec.mua).^2);
   metric_eval.rmse = sqrt(MSE);
   peakVal = max(mesh_rec.mua);
   metric_eval.psnr = 10*log10(peakVal*peakVal/MSE); 
   signal_front = sum(mesh_rec.mua.^2);
   noise_back = sum((mesh_GT.mua - mesh_rec.mua).^2);
   metric_eval.snr = 10*log10(signal_front./noise_back);
   
    % segmentation
    meshSeg_GT = mesh_GT;
    mesh_segmented = mesh_rec;
    if isnan(thresh)
        thresh = 0.5*(max(mesh_rec.mua) + median(mesh_rec.mua));       
        thresh_GT = 0.5*(max(mesh_GT.mua) + median(mesh_GT.mua));
    
%         thresh =0.5*peakVal;  
%         thresh_GT = 0.5*max(mesh_GT.mua);
        
        mesh_segmented.mua = double(mesh_rec.mua > thresh);
        meshSeg_GT.mua = double(mesh_GT.mua > thresh_GT);
    end
    % evaluate the structual reconstruction
    metric_eval.dice = dice(...
           mesh_segmented.mua, meshSeg_GT.mua);
    
     % evaluate resolution (to be coded)
    if isResolution
       metric_eval.resolution = 0;
    end
else
    disp( 'evaluate mus to be coded')
end