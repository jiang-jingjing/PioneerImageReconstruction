%>@brief Brief description of the function
%>
%> prepare data from measurement for time domain modality 
%>
%>@param  d diameter
%>@param  z depth
%>@param  unitmm, unit in mm 
%>
%> @retval volOut create volume
  
function [volOut, varargout] = createCylinderMC(d, varargin)
% default
is2D = 1;
unitmm = 1; 
isSave = 0;
if ~isempty(varargin)
    if length(varargin) >= 1
        % log file name
        if numel(varargin{1})  
            if  varargin{1} > 0
                z = varargin{1};
                is2D = 0;
            else
                warning(['unexpected non-positive number for dimension z']);
            end
        else
            error('Bad 2th argument value. number expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 2
        % if differential imaging
        if numel(varargin{2}) 
            if  varargin{2} > 0
                unitmm = varargin{2};
            else
                warning(['unexpected non-positive number for unitmm']);
            end

        else
            error('Bad 3th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
     if length(varargin) >= 3
        % save volume to .mat
        if ischar(varargin{3}) 
            
                file_vol = varargin{3};
                isSave = 1;
        else
            error('Bad 4th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
end

 
    
if is2D
    [X Y ] = meshgrid(1:d);
    R = d/2;
    mask_cy = (X-R-0.5).^2 + (Y-R-0.5).^2 < (R-1).^2 ;
    volOut =  uint8(mask_cy);
     
    warning('2D MC is not tested')
else  
    [X Y Z] = meshgrid(1:d, 1:d,1:z);
    R = d/2;
    mask_cy = (X-R-0.5).^2 + (Y-R-0.5).^2 < (R-1).^2 ;
%     volOut = zeros(d,d,z);
    volOut =  uint8(mask_cy);
 end

if isSave
    save(file_vol, 'volOut')
end