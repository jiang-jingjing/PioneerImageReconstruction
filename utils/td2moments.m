%>@brief Brief description of the function
%>
%> convert time domain data to moments
%>
%>@param dataFwd forward result
%>          .time
%>          .tpsf

%>@param n_moments  maximum order + 1
%>@param paras parameters  
%>@param cfg  struct config parameter for time gates 
%>          .tstart
%>          .tend
%>          .tstep 
%>          .prop optical properties
%> 
%> @retval dataFwd forward result updated
%> author: jingjing jiang jing.jing.jiang@outlook.com
function dataFwd = td2moments(dataFwd, n_moments)
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

dataFwd.moments = moments;
dataFwd.n_moments = n_moments;
end


 