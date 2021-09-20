%>@brief Brief description of the function
%>
%> image reconstruction for time domain modality 
%>
%>@param volIn tissue volume
%>@param paras parameters  
%>
%> @retval volOut reconstructed volume
%> author: jingjing jiang jing.jing.jiang@outlook.com
function [volOut, varargout] = reconstructionTime(volIn, paras, varargin)
    
for iter = 1:num_iter
    %% calculate jacobian 
    
    
    
    
    %% calculate error
    data_diff = calculate_data_difference(value_rec, data_meas);

    %% reconstruction
    d=size(J);
    L=eye(d(1,2),d(1,2));
    tol = 1e-8;
    maxit = 1e3;

    % 
    [value_rec,fl,rr,it,rv,lsrv]=lsqr([...
        J;reg*L],[data_diff;zeros(d(1,2),1)],...
        tol,maxit,[],[],[]);
    it 

    %% add the changes to the volume
    vol_rec = update_vol(vol_rec, value_rec, paras);
    
    %% stopping criteria
    
end
end

function [data_diff] = calculate_data_difference(data_update, data_measured)
    num_gates = size(data_measured.tpsf,2);
    reference =   log(data_update.tpsf +...
                    (data_update.tpsf <=0));
    anomaly = log(data_measured.tpsf+...
                (data_measured.tpsf<=0));
    if num_gates > 1
        for ior = 2:num_gates
            reference = [reference; log(data_update.tpsf(mask_data, ior)+...
        (data_update.tpsf(mask_data, ior)<=0));];
             anomaly = [anomaly; log(data_measured.tpsf(mask_data, ior)+...
        (data_measured.tpsf(mask_data, ior)<=0))];
        end
    end
            
    % split to match the Jacobian interlaced format
    len = size(reference,1);
    data_diff = zeros(len,1);
    mm = len/num_gates;
    for ior = 1 : num_gates
        data_diff(ior:num_gates:len) = ...
            reference(1+mm*(ior-1):mm*ior) ...
            - anomaly(1+mm*(ior-1):mm*ior);
    end
end