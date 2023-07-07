%>@brief Brief description of the function
%>
%> calculate instrumental response function for TD systems
%> getIRF.m
%> [gt g_fft] = getIRF(fMeas, fSim, time)
%> @param fMeas: tpsf from measurement
%> @param fSim: tpsf from simulation
%> @param time: time s
%> @param freqT [optional]: remove high frequencies>freqT
%>  
%> @retval gt IRF in time domain
%> @retval g_fft IRF in frequency domain 
%> author: jingjing jiang jing.jing.jiang@outlook.com

 
function [gt g_fft] = getIRF_TD(fMeas, fSim, time,varargin)
    freqT =0;
    % varagin
    if ~isempty(varargin)
       if length(varargin) >= 1
               % interpolation :flagInterp= 1
            if isnumeric(varargin{1}) 
                freqT = varargin{1};
            else
                error('Bad 4th argument value. number expected for freqT.')
            end  
       end
    end
    numTG = length(time);

    N = 2048;
    tstep = median(diff(time));
    f = (1:N/2)  / N / tstep; % tstep = 0.1 ns

    y_fft_ref =fft(fMeas, N) ;
    y_nirfast_fft = fft(fSim,N);
    g_fft = y_fft_ref./y_nirfast_fft;
     
    if freqT>0 % low pass
        id_th = find(f<freqT);
        iis = id_th(end);
        g_fft([iis+1:N/2 (N/2+1):(N-iis)]) = 0;
    end
    figure(101), subplot(121)
    plot(abs(g_fft))
    title('FT IRF')
    gt = ifft(g_fft, N);
    
    gt = gt(1:numTG);
    figure(101), subplot(122)
    plot(time,abs(gt))
    title('IRF')
end
