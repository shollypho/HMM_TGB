%% renaming script to save logfiles with movement information
cd /imaging/hp02/TGB/rest_closed/MF/Controls/
con_data = Control_Subs_Data;
for ss = 1:length(con_data)
    ss
    
    folder_s1 = con_data{1,ss}.Name{1};
    folder_s2 = con_data{1,ss}.Name{2};
    
    
    
    if ~strcmp(folder_s1, 'NA')
        cd(folder_s1)
        sss_file_s1 = sprintf('%s_1_ssslogfile.log', folder_s1);
        if exist(sss_file_s1, 'file')
            copyfile(sss_file_s1, sprintf('movecomp_%s', sss_file_s1));
        end
        
        sss_file_s2 = sprintf('%s_2_ssslogfile.log', folder_s1);
        
        if exist(sss_file_s2, 'file')
            copyfile(sss_file_s2, sprintf('movecomp_%s', sss_file_s2));
        end
        cd ..
    end
    
    
    
    if ~strcmp(folder_s2, 'NA')
        cd(folder_s2)
        
        sss_file_s1 = sprintf('%s_1_ssslogfile.log', folder_s2);
        if exist(sss_file_s1, 'file')
            copyfile(sss_file_s1, sprintf('movecomp_%s', sss_file_s1));
        end
        
        
        sss_file_s2 = sprintf('%s_2_ssslogfile.log', folder_s2);
        
        if exist(sss_file_s2, 'file')
            copyfile(sss_file_s2, sprintf('movecomp_%s', sss_file_s2));
        end
        cd ..
    end
    
    
end