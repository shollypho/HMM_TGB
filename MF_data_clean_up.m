%% MF clean up
% load control and patient data
con_data = Control_Subs_Data;
pat_data = Patient_Subs_Data;

% path to all maxfiltered files
mf_path = '/imaging/hp02/TGB/rest_closed/MF';
cd(mf_path)  

% remove excess control data - only need the rsc#_trans_meg###.fif s
cd 'Controls'
con_trans_check = [];
for ss = 1:length(con_data) % subjects
    ss
    for k = 1:2 % sessions
        % check if session exists (important for patients)
        if ~strcmp(con_data{1,ss}.Name{k}, 'NA')
            cd(con_data{1,ss}.Name{k})
            
            for blk = 1:2 % blocks
                % check if trans file exists
                if exist(sprintf('rsc%d_trans_%s.fif',blk,con_data{1,ss}.Name{k}), 'file')
                    
                    delete(sprintf('rsc%d_sss_%s.fif',blk,con_data{1,ss}.Name{k}))
                    delete(sprintf('rsc%d_sss_%s_wods.fif',blk,con_data{1,ss}.Name{k}))
                    delete(sprintf('%s_rsc%d_badCH.fif',con_data{1,ss}.Name{k},blk))
                    
                    if blk == 2
                        delete(sprintf('rsc%d_transA_%s.fif',blk,con_data{1,ss}.Name{k}))
                    end
                    
                    con_trans_check((ss*blk),k) = 1;
                else
                    con_trans_check((ss*blk),k) = 0;
                end
                
                
            end

        end
        cd(mf_path)
        cd 'Controls/'
    end
end



%% Now for patients:

% path to all maxfiltered files
mf_path = '/imaging/hp02/TGB/rest_closed/MF';
cd(mf_path)  

% remove excess control data - only need the rsc#_trans_meg###.fif s
cd 'Patients'
pat_trans_check = [];
for ss = 1:length(pat_data) % subjects
    ss
    for k = 1:2 % sessions
        % check if session exists (important for patients)
        if ~strcmp(pat_data{1,ss}.Name{k}, 'NA')
            cd(pat_data{1,ss}.Name{k})
            
            for blk = 1:2 % blocks
                % check if trans file exists
                if exist(sprintf('rsc%d_trans_%s.fif',blk,pat_data{1,ss}.Name{k}), 'file')
                    
                    delete(sprintf('rsc%d_sss_%s.fif',blk,pat_data{1,ss}.Name{k}))
                    delete(sprintf('rsc%d_sss_%s_wods.fif',blk,pat_data{1,ss}.Name{k}))
                    delete(sprintf('%s_rsc%d_badCH.fif',pat_data{1,ss}.Name{k},blk))
                    
                    if blk == 2
                        delete(sprintf('rsc%d_transA_%s.fif',blk,pat_data{1,ss}.Name{k}))
                    end
                    
                    pat_trans_check((ss*blk),k) = 1;
                else
                    pat_trans_check((ss*blk),k) = 0;
                end
                
                
            end

        end
        cd(mf_path)
        cd 'Patients/'
    end
end