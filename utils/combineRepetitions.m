function combineRepetitions(fldr,waveList, group_list, ...
    rep_list, srcList,Prefix, isPM, scalerPM)
 
for i_wav = 1:length(waveList)
    wav = num2str(waveList(i_wav))
    for i_src = srcList
        for i_group = 1:length(group_list)
            group = group_list{i_group};
            for rep_no = rep_list
                foo = sprintf([fldr '/' 'timing_data_corrected/' group '/timing_response_%s_%s_%s_%s.mat'], ...
                    Prefix, wav, num2str(i_src), num2str(rep_no));
                load(foo, 'timing_response')
                if rep_no==1
                    if isPM
                        f_pm= sprintf([fldr '/' 'timing_data_corrected/' group '/powermeter_%s_%s_%s_%s.txt'], ...
                        Prefix, wav, num2str(i_src), num2str(rep_no));
                        formatSpec = '%f';
                        fileID = fopen(f_pm, 'r');
                        pw =  fscanf(fileID, formatSpec);
                        pw = nanmean(pw) .* scalerPM;
                        
                        timing_response = ...
                            double(timing_response)./pw;         
                    end
                    data(i_group).timing_response = double(timing_response);
                else
                    if isPM
                        f_pm= sprintf([fldr '/' 'timing_data_corrected/' group '/powermeter_%s_%s_%s_%s.txt'], ...
                        Prefix, wav, num2str(i_src), num2str(rep_no));
                        formatSpec = '%f';
                        fileID = fopen(f_pm, 'r');
                        try
                            pw =  fscanf(fileID, formatSpec);
                        catch
                            disp(['no pw value at repetition ' ...
                                num2str(rep_no) ', use the previous repetition'])
                            
                        end
                        pw = nanmean(pw) .* scalerPM;
                        
                        timing_response = ...
                            double(timing_response)./pw;         
                    end
                    data(i_group).timing_response = ...
                        data(i_group).timing_response +...
                        double(timing_response);
                    
                end
                
            end
           timing_response = data(i_group).timing_response;
           foo_allrep = sprintf([fldr '/' 'timing_data_corrected/' group '/timing_response_%s_%s_%s.mat'], ...
                    Prefix, wav, num2str(i_src));
            save(foo_allrep, 'timing_response') 
        end
    end
end