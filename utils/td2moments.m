%>@brief Brief description of the function
%>
%> convert time domain data to moments
%>
%>@param dataFwd forward result
%>          .time
%>          .tpsf

%>@param n_moments  maximum order + 1
%>@param varargin: calibration
%>              @param dataMEAS_ref
%>              @param dataSIM_ref
%> 
%> @retval dataFwd forward result updated
%> author: jingjing jiang jing.jing.jiang@outlook.com
function dataFwd = td2moments(dataFwd, n_moments, varargin)
isExistRef = 0;
if ~isempty(varargin)
   if length(varargin) >= 1
       if isstruct(varargin{1}) 
            dataMEAS_ref = varargin{1};
            isExistRef = 1;
       else
            error('Bad 3th argument value. struct expected for dataMEAS_ref.')
       end
   end
    if length(varargin) >= 2
        if isstruct(varargin{2}) 
            dataSIM_ref = varargin{2};
        else
            error('Bad 4th argument value. struct expected for dataSIM_ref.')
        end
    end 
end
tstep = unique(diff(dataFwd.time));
if length(tstep)>1
    if std(tstep)>1e-12
    warning('time step is not the same') 
    end
%     tstep = mean(tstep);    
end
time = dataFwd.time; 
tof_list = dataFwd.tpsf;
count = size(tof_list,1);
moments = zeros(count, n_moments);
for ind = 1:count
    tof = tof_list(ind,:);
    mm = zeros(n_moments,1);
    mm(1) = sum(tof);
    if n_moments>1
         mm(2) = (time*tof')/mm(1);
    end
    if n_moments > 2
        for ii = 2:(n_moments-1)
            mm(ii + 1) = (tof*(time').^ii)/mm(1);
        end
    end
    moments(ind, :) = mm;
end 

%% calibration
if isExistRef
    im = 1;
    scaler  = dataMEAS_ref.moments(:,im) ./ ...
    dataSIM_ref.moments(:,im);
    data_cal.moments(:,im) = moments(:,im)./ ...
        scaler; 
    for im = 2:size(moments,2)
        scaler  = dataMEAS_ref.moments(:,im) - ...
            dataSIM_ref.moments(:,im);
        data_cal.moments(:,im) = moments(:,im) - ...
            scaler;
    end
    dataFwd.moments = data_cal.moments;
    dataFwd.link = dataSIM_ref.link;
else
    
dataFwd.moments = moments;

end
dataFwd.n_moments = n_moments;
end


 