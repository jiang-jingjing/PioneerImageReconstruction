%>@brief Brief description of the function
%>
%> convert time domain data to Fourier domain
%>
%>@param dataFwd forward result
%>          .time
%>          .tpsf

%>@param freq  frequency (Hz)  
%>@param paras parameters  
%>@param cfg  struct config parameter for time gates 
%>          .tstart
%>          .tend
%>          .tstep 
%>          .prop optical properties
%> 
%> @retval dataFwd forward result updated
%> author: jingjing jiang jing.jing.jiang@outlook.com
function dataFwd = td2fd(dataFwd, freq)
tstep = unique(diff(dataFwd.time));
if length(tstep)>1
    if std(tstep)>1e-12
    warning('time step is not the same') 
    end
    tstep = mean(tstep);    
end
N = 2048;
f = (1:N/2)' / N / tstep; % tstep = 12.5ns/256
numfreq = length(freq);

count = size(dataFwd.tpsf,1); %count: number of histograms
tof_list = dataFwd.tpsf;
for iff = 1:numfreq
%     df = dfList(iff);%  
    df = find(abs(freq(iff) - f)<5e6);
    for tt = 1:  count  
         temp =fft(tof_list(tt,:), N) ;
         temp_fft = temp(df);
         amp = abs(temp_fft);
         amp_log =   log(temp_fft );
         ph_temp = abs(unwrap(angle(temp_fft)));
         data_fft.amplitude(tt) = amp; % list of amplitude
         data_fft.phase(tt) =  ph_temp.*180 /pi; % list of phase
    end 
    
    if iff==1
        data_fft_all = data_fft;
    else
        data_fft_all.amplitude = [data_fft_all.amplitude ...
            data_fft.amplitude];
        data_fft_all.phase = [data_fft_all.phase ...
            data_fft.phase];
    end
end

dataFwd.amplitude = data_fft_all.amplitude';
dataFwd.phase = data_fft_all.phase';
dataFwd.frequency = freq;
end
