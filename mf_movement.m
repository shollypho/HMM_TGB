% Record Translation and Rotation movements to assess the amount of 
% movement throughout the resting block

con_data = Control_Subs_Data;
pat_data = Patient_Subs_Data;

% path to all maxfiltered files
mf_path = '/imaging/hp02/TGB/rest_closed/MF';
cd(mf_path)  
cd Patients
result = [];
j = 1;
for ss = 5:length(pat_data) % subjects
    ss
    for k = 1:2 % sessions
        % check if session exists (important for patients)
        if ~strcmp(pat_data{1,ss}.Name{k}, 'NA')
            cd(pat_data{1,ss}.Name{k})
            
            for blk = 1:2 % blocks
                
                %if ss == 17 && blk == 2 % problem with this block in control 17
                %    continue
                % check if trans file exists
                if exist(sprintf('movecomp_%s_%d_ssslogfile.log',pat_data{1,ss}.Name{k},blk), 'file')
                    
                    [linet, lineg, linev, liner, lined] = record_movecomp(sprintf('movecomp_%s_%d_ssslogfile.log',pat_data{1,ss}.Name{k},blk));
                    
                    %Search for time HPIs are on:
                    move{ss}.hpi_time(k,blk) = sum(lineg>0.99);
                    % times translation is more than 0.3cm/s:
                    move{ss}.trans_move(k,blk) = sum(linev>0.3);
                    % times rotation is more than 0.3rads/s:
                    move{ss}.rotat_move(k,blk) = sum(liner>0.3);
                    % times drift is more than 3mm:
                    move{ss}.drift_move(k,blk) = sum(lined>0.3);
                    
                    %remaining time:
                    move{ss}.good_move(k,blk) = move{ss}.hpi_time(k,blk)-move{ss}.trans_move(k,blk)-move{ss}.rotat_move(k,blk);
                   
                    result{j,1} = pat_data{1,ss}.Name{k};
                    result{j,2} = k;
                    result{j,3} = blk;
                    result{j,4} = move{ss}.hpi_time(k,blk);
                    result{j,5} = move{ss}.trans_move(k,blk);
                    result{j,6} = move{ss}.rotat_move(k,blk);
                    result{j,7} = move{ss}.drift_move(k,blk);
                    result{j,8} = move{ss}.good_move(k,blk);
                    
                    good_move_summ{k}{ss,1} = pat_data{1,ss}.Name{k};
                    good_move_summ{k}{ss,blk+1} = move{ss}.good_move(k,blk);
                    
                    j= j+1;
                end
            end
            cd(mf_path)  
            cd Patients
        end
    end
end

