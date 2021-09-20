%>@brief Brief description of the function
%>
%> calculate forward result for time domain modality 
%> based on monte carlo simulation
%>
%>@param volIn tissue volume
%>@param paras parameters  
%> 
%> @retval h handle for figure
%> author: jingjing jiang jing.jing.jiang@outlook.com

function   h =plotVolMC(volIn, varargin)
plotSrcDet = 0;
if ~isempty(varargin)
    if length(varargin) >= 1
        % log file name
        if isstruct(varargin{1})  
          paras = varargin{1};
                plotSrcDet = 1;
            
           
        else
            error('Bad 2th argument value. struct expected. Please see help for details on how to use this function.')
        end
    end
end
   [  h] = mcxplotvol(volIn);
    

    colorbar;
    if plotSrcDet 
        det = paras.det.coord;
        src = paras.src.coord;
        hold on
        plot3(det(:,1), det(:,2), ...
            det(:,3), 'ob' )
        plot3(src(:,1), src(:,2), ...
            src(:,3), 'or' )
        plot3(src(1,1), src(1,2), ...
            src(1,3), 'og' )
%         legend('detetors','sources','source 1')
    end
        
end