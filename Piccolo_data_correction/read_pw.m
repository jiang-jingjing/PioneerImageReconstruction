pw_all = [];
tPW_all = [];
timepath = ['/Timestamps/' ];
formatIn = 'HH:MM:SS';
time_meas_str = h5read(file, timepath) ;
time_meas_num = [datenum(time_meas_str(:,1),formatIn), ...
    datenum(time_meas_str(:,2),formatIn)] .* 24*3600;

 
for i_id = 1:length(ids)
    pwpath = ['/Power_log/' num2str(ids(i_id))];
    pwTpath = ['/Power_times/'  num2str(ids(i_id))];
try    pw_value = double(h5read(file,pwpath));
catch
    disp(['no pw value at ' num2str(i_id) 'th data point'])
    continue
end
pw_time = double( h5read(file, pwTpath)) ;
    
    pw_all = [pw_all; pw_value];
    tPW_all = [tPW_all; pw_time./1e3  + time_meas_num(ids(i_id),1)];
    
    
    datapath = ['/Data/' num2str(ids(i_id))];
    wav = num2str(h5readatt(file, datapath, 'Wavelength [nm]'));
    source_no = h5readatt(file, datapath, 'Channel');
    rep_no = h5readatt(file, datapath, 'Run');
    target_pos = h5readatt(file, datapath, 'Targetted Head position {x, y, z, dist}');
   
 
    
    file_pos = ['target_' num2str(target_pos(1)) '_' ...
        num2str(target_pos(2)) '_' ...
        num2str(target_pos(3)) '_' ...
        num2str(target_pos(4)) ];
    if exist([fldr '/timing_data_corrected/' file_pos])==0
        mkdir([fldr '/timing_data_corrected/' file_pos])
    end
%     
%     
fout = sprintf([fldr '/' 'timing_data_corrected/' file_pos '/powermeter_%s_%s_%s_%s.txt'], ...
    Prefix, wav, num2str(source_no), num2str(rep_no));
fid = fopen(fout, 'w');
fprintf(fid, '%.10f\n', pw_value);
fclose(fid);
end

figure,subplot(121)
plot(tPW_all-min(tPW_all)  , pw_all,'.-')
title(['powermeter ' num2str(waveList) 'nm'])
xlabel('seconds')
subplot(122), boxplot(pw_all)