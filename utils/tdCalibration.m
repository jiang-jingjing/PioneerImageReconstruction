%>@brief Brief description of the function
%>
%> remove instrumental response function
%> IRF is calculated with a measurement on homogeneous phantom
%>
%>@param dataFwd forward result
%>          .time
%>          .tpsf

%>@param time:  time gates  
%>@param dataMEAS_ref
%>@param dataSIM_ref
%> varargin: 
%>
%> @retval dataFwd forward result updated
%> author: jingjing jiang jing.jing.jiang@outlook.com
function [data_fwd] = tdCalibration(dataFwd, ...
                time, dataMEAS_ref, dataSIM_ref,varargin)

flagInterp = 0;
SmoothDelta = 0;
% varagin
if ~isempty(varargin)
   if length(varargin) >= 1
           % interpolation :flagInterp= 1
        if isnumeric(varargin{1}) 
            flagInterp = varargin{1};
        else
            error('Bad 5th argument value. number expected for flagInterp.')
        end  
   end
   if length(varargin) >= 2
           % smoothing
        if isnumeric(varargin{2}) 
            SmoothDelta = varargin{2};
        else
            error('Bad 6th argument value. number expected for SmoothDelta.')
        end  
   end
end
 
% The relative tolerance stays the same.
OPTIONS = solver_options;
OPTIONS.tolerance = 1e-12;
%% prepared NIRFAST forward result for meaured data
tof_list = dataFwd.tpsf;
tof_list_ref = dataMEAS_ref.tpsf;
 
t_nirfast= dataSIM_ref.time;


% -------- interpolation x2 num gates
if flagInterp % time resolution: 0.0488/2 instead of 0.0488
[rr cc] = size(tof_list);
tof_list_2x = zeros(rr, cc*2);
tof_list_ref_2x = tof_list_2x;
 for id = 1:rr
tof_list_2x(id,:) = interp1(time, tof_list(id, :), t_nirfast, 'linear');
tof_list_ref_2x(id,:) = interp1(time, tof_list_ref(id, :), t_nirfast, 'linear');

% figure, 
% semilogy(time,tof_list(id,:),'o-')
% hold on
% plot(t_nirfast, tof_list_2x(id,:),'x-')
% semilogy(time,tof_list_ref(id,:),'o-')
% plot(t_nirfast, tof_list_ref_2x(id,:),'x-')
 end

tof_list = tof_list_2x;
tof_list_ref = tof_list_ref_2x;
tof_list( isnan(tof_list)) =  0;
tof_list_ref( isnan(tof_list_ref)) =  0;

time = t_nirfast;
clear tof_list_2x tof_list_ref_2x
%------------------ interpolation end
end

%% convolution
%% FFT
numTG = length(time);

N = 2048;
tstep = median(diff(time));
f = (1:N/2)' / N / tstep; % tstep = 12.5ns/256
count = size(tof_list,1); %count: number of histograms

tof_cal_list = zeros(size(tof_list));
for tt = 1:  count  

%     % filtering
% d1 = designfilt('lowpassiir','FilterOrder',12, ...
%     'HalfPowerFrequency',0.15,'DesignMethod','butter');
% x = tof_list_ref(tt,:);
% tof_ref_filt = filtfilt(d1,x);
% 
% figure,
% semilogy(x)
% hold on
% semilogy(tof_ref_filt)
% y_fft_ref =fft(tof_ref_filt, N) ;

     y_fft_ref =fft(tof_list_ref(tt,:), N) ;
     y_nirfast_fft = fft(dataSIM_ref.tpsf(tt,:),N);
     g_fft = y_fft_ref./y_nirfast_fft;
 
     y_fft =fft(tof_list(tt,:), N) ;
     y_fft_cal = y_fft ./ g_fft;
      
     temp_cal = ifft(y_fft_cal, N);
     tof_cal =temp_cal(1:numTG);
%      figure,
%      subplot(121)
%      semilogy(tof_list_ref(tt,:))
%      hold on
%      plot(tof_list(tt,:))
%      legend('meas homo', 'meas heter')
%      subplot(122)
%      semilogy(dataSIM_ref.tpsf(tt,:))
%      hold on
%    
%      plot(tof_cal)
%      legend('nirfast homo', 'meas heter calibrated') 
     
     tof_cal_list(tt,:) = tof_cal;
end 
 data_fwd = dataSIM_ref;
 data_fwd.tpsf = tof_cal_list;

end


function dataFilt = lowpassIIR( data)
d = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.4,'DesignMethod','butter');
dataFilt = filtfilt(d, data);
  
figure(155),
 
cla
semilogy( data)
hold on
plot( dataFilt)
 
legend('raw','filtered')
title('lowpass filtered')
end

function dataFilt = myfilter(SmoothDelta, data)

% sigma = 1;
% sz = 50;    % length of gaussFilter vector
% x = linspace(-sz / 2, sz / 2, sz);
% gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
% gaussFilter = gaussFilter / sum (gaussFilter); % normalize

%     dataFilt.amplitude = filter(gaussFilter,SmoothDelta,...
%         data.amplitude);
%     dataFilt.phase = filter(gaussFilter,SmoothDelta,...
%         data.phase);
 
dataFilt.amplitude = medfilt1(data.amplitude,SmoothDelta);
dataFilt.phase = medfilt1(data.phase,SmoothDelta);
 
figure(155),
 
subplot(121)
cla
semilogy(data.amplitude)
hold on
plot(dataFilt.amplitude)
subplot(122)
cla
plot(data.phase)
hold on
plot(dataFilt.phase)
legend('raw','filtered')
end