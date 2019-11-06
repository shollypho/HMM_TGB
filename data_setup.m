function all_data = data_setup()

% Load Control and Patient data names
Con_Sub = Control_Subs_Data;
Pat_Sub = Patient_Subs_Data;

load('/imaging/hp02/TGB/rest_closed/MF/hpi_turned_on.mat');
%mf_HPIs_on;

PSP_match_Cons = {'C1', 'C2','C3','C4','C8','C11', 'C12','C13','C14','C15'};
bv_match_Cons = {'C5', 'C6','C7','C9','C16','C17', 'C18','C19','C20'};



j = 1;
all_data = {};
for ss = 1:length(Con_Sub)
    if sum(Con_Sub{ss}.RSC) ~= 0 % if the participant has at least one good run
        for k = 1:4 % loop through the 2 sessions and 2 runs
            if Con_Sub{ss}.RSC(k)==1
                all_data{j,1} = Con_Sub{ss}.Name{3}; % ID
                all_data{j,2} = 'Con'; % Control group
                all_data{j,3} = Con_Sub{ss}.Name{ceil(k/2)}; % MEG ID
                all_data{j,4} = ceil(k/2); % Session
                all_data{j,5} = k; % run
                if k ==3; all_data{j,5}=1; elseif k == 4; all_data{j,5}=2; end
                all_data{j,6} = Con_Sub{ss}.Name{4}; % MRI ID
                all_data{j,7} = Con_Sub{ss}.DrugSess==ceil(k/2); % Drug (1), Placebo = (0)
                
                % extract fids
                if ~isempty(Con_Sub{ss}.Fids)
                    Fids = Con_Sub{ss}.Fids';
                    for f = 1:9
                        all_data{j,7+f} = Fids(f);
                    end
                end
                
                % fif file path, name, exist
                all_data{j,17} =['/imaging/hp02/TGB/rest_closed/MF/Controls/',all_data{j,3}, '/' ];
                all_data{j,18} =sprintf('rsc%d_sss_%s.fif', all_data{j,5}, all_data{j,3});
                %all_data{j,19} = exist([all_data{j,17},all_data{j,18}],'file');
                
                % New spm file name
                all_data{j,19} = sprintf('%s_s%d_rsc%d',all_data{j,1}, all_data{j,4}, all_data{j,5});
                
                % EEG recorded:
                all_data{j,20} = Con_Sub{ss}.EEG(all_data{j,4});
                
                % HPI turned on in 250Hz
                all_data{j,21} = hpi_start_250hz(j);
                
                 
                    for cc = 1:length(PSP_match_Cons)
                        if strcmp(all_data{j,1},PSP_match_Cons{cc})
                           all_data{j,23} = 2; 
                        end
                    end
                     for cc = 1:length(bv_match_Cons)
                        if strcmp(all_data{j,1},bv_match_Cons{cc})
                           all_data{j,23} = 4; 
                        end
                    end
                
                j= j+1;
            end
        end
    end
    
end

for ss = 1:length(Pat_Sub)
    if sum(Pat_Sub{ss}.RSC) ~= 0 % if the participant has at least one good run
        for k = 1:4 % loop through the 2 sessions and 2 runs
            if Pat_Sub{ss}.RSC(k)==1
                all_data{j,1} = Pat_Sub{ss}.Name{3}; % ID
                all_data{j,2} = Pat_Sub{ss}.Dia; % Patient groups
                all_data{j,3} = Pat_Sub{ss}.Name{ceil(k/2)}; % MEG ID
                all_data{j,4} = ceil(k/2); % Session
                all_data{j,5} = k; % run
                if k ==3; all_data{j,5}=1; elseif k == 4; all_data{j,5}=2; end
                all_data{j,6} = Pat_Sub{ss}.Name{4}; % MRI ID
                all_data{j,7} = Pat_Sub{ss}.DrugSess==ceil(k/2); % Drug (1), Placebo = (0)
                
                % extract fids
                if ~isempty(Pat_Sub{ss}.Fids)
                    Fids = Pat_Sub{ss}.Fids';
                    for f = 1:9
                        all_data{j,7+f} = Fids(f);
                    end
                end
                
                % fif file path, name, exist
                all_data{j,17} =['/imaging/hp02/TGB/rest_closed/MF/Patients/',all_data{j,3}, '/' ];
                all_data{j,18} =sprintf('rsc%d_sss_%s.fif', all_data{j,5}, all_data{j,3});
                %all_data{j,19} = exist([all_data{j,17},all_data{j,18}],'file');
                
                % New spm file name
                all_data{j,19} = sprintf('%s_s%d_rsc%d',all_data{j,1}, all_data{j,4}, all_data{j,5});
                
                % EEG recorded:
                all_data{j,20} = Pat_Sub{ss}.EEG(all_data{j,4});
                
                % HPI turned on in 250Hz
                all_data{j,21} = hpi_start_250hz(j);
                
                
                
                % Separate into patients and matched control groups
                % PSP = 1, Control matched with PSP = 2, bvFTD = 3, Control
                % matched with bvFTD = 4
                if strcmp(all_data{j,2},'PSP')
                    all_data{j,23} = 1;
                elseif strcmp(all_data{j,2},'bvFTD')
                    all_data{j,23} = 3;
               
                end
                
                
                
                j= j+1;
            end
        end
    end
end

%% Add in which participant scans are now included
for j = 1:size(all_data,1)
    % Exclude second rsc
    if all_data{j,5} == 2 
        all_data{j,22} = 0;
        % Exclude outliers
    elseif  j==81 || j==82 || j==35 || j==37 || j==119 || j == 120 || j == 131 || j==132%||j==95 || j==100 || j==131 || j==14 
        all_data{j,22} = 0;
        % Exclude any without both sessions
    elseif j == 152 || j == 144 ||j==130 || j==96 || j==98 %j==13 || j==82 || j==94|| j==101 || j==120 || j==132 
        all_data{j,22} = 0;
    elseif j == 94 || j == 95 || j == 100 || j == 101% P13 & p17, too short after movement removal
        all_data{j,22} = 0;
    else
        all_data{j,22} = 1;
    end
end