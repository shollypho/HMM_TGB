%% GLEAN pipeline

%% Setup OSL

addpath(genpath('/imaging/hp02/software_n_scripts/OSL/osl/osl-core'))
osl_startup
osl_check_installation
open osl_debug_log.txt

%% Define paths

addpath('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/scripts')

data_dir = '/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/CamCAN_F_provisional_analysis_fids_and_headshape_no_EEG_128_no_template'; % path where data is located
glean_name = '/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/CamCAN_F_provisional_analysis_fids_and_headshape_no_EEG_128_no_template/GLEAN/GLEAN.mat'; % name for this GLEAN analysis
glean_dir = fullfile(osldir,'example_data','glean_example'); % path for GLEAN stuff (masks etc) 
check_coreg_directory_EEG = '/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/CamCAN_F_provisional_analysis_fids_and_headshape_no_EEG_128_no_template/Check_coreg_fiducials_figures_output_EEG/';
check_coreg_directory_MEG = '/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/CamCAN_F_provisional_analysis_fids_and_headshape_no_EEG_128_no_template/Check_coreg_fiducials_figures_output_MEG/';

cd(data_dir)

nStates = 8; % number of states to be used for GLEAN
nStates_str = num2str(nStates);

%% Get data

q = {};

[NUMERIC TXT D] = xlsread('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/CamCAN_F_MRI_Series_Nos_with_IDs_updated_for_Roni_HMM.xlsx'); % read data from xls file. David - you can add the full path, or CD to the directory where the file is located
megPath = '/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/Maxfiltered_all_MEG_tasks_and_resting_state/';
mriPath = '/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/aamod_coreg_extended_1_00001_with_fixed_CC710758_and_7_needing_template/';
fidPath = '/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/';

q.SubCCIDc = D(2:end,1)'; % D(2:end,1)% get CCIDs for all subjects
megFolder = D(2:end,2)'; % get MEG folder for all subjects

for s = 1:length(q.SubCCIDc)
    
    if ~strcmp(q.SubCCIDc{s},'CC510371') && ~strcmp(q.SubCCIDc{s},'CC620001') && ~strcmp(q.SubCCIDc{s},'CC710428') && ~strcmp(q.SubCCIDc{s},'CC720167') && ~strcmp(q.SubCCIDc{s},'CC721184') && ~strcmp(q.SubCCIDc{s},'PP138368') && ~strcmp(q.SubCCIDc{s},'PP185442') ==1;    
        
        %% Get MEG data
    
        try
            q.FileNames.meg{1,s} = strtrim(ls(sprintf('%s%s/rest/tsss.fif',megPath,lower(megFolder{s}))));
            chkFile(1,s) = 1;
        catch
            warning('No MEG file for subject %s',q.SubCCIDc{s}); 
            q.FileNames.meg{1,s} = [];
        end

        %% Get MRI data

        try
            q.FileNames.mri{1,s} = strtrim(ls(sprintf('%s%s/structurals/*.nii',mriPath,q.SubCCIDc{s})));
            chkFile(2,s) = 1;
        catch
            warning('No MRI file in mriPath for subject %s',q.SubCCIDc{s}); 
            q.FileNames.mri{1,s} = [];
        end


        %% Get Fiducials data

        try
            q.FileNames.fid{1,s} = strtrim(ls(sprintf('%s%s/FiducialLocs.mat',fidPath,q.SubCCIDc{s})));
            chkFile(3,s) = 1;
        catch
            warning('No FID file for subject %s',q.SubCCIDc{s}); 
            q.FileNames.fid{1,s} = [];
        end
    end
end;

q.AllExistSubInd = find(sum(chkFile)==3); % get indices for subjects that have all files: MEG, MRI, FID

save(sprintf('%s/q',data_dir),'q')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coregistration and Forward model

load(sprintf('%s/q',data_dir))

fileindex = q.AllExistSubInd; % those subjects with both anatomical and fiducial markers
subs = q.SubCCIDc(fileindex);
subjects = subs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start MATLAB pool
openparallel = 1;


 if openparallel && numel(gcp('nocreate')) == 0
     %if matlabpool('size')==0
         MaxNsubs = max([16 length(subjects)]);
         P = cbupool(MaxNsubs);
         P.ResourceTemplate = '-l nodes=^N^,mem=4GB,walltime=140:00:00';
         P.SubmitArguments='-W x="NODESET:ONEOF:FEATURES:MAXFILTER"';
         poolhandle = parpool(P);
     %end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

cd(sprintf('%s/MEG',data_dir))

parfor s = 1:length(subs)
     
    destMEG = sprintf('%s/MEG/%s.fif',data_dir,q.SubCCIDc{fileindex(s)});
    destMRI = sprintf('%s/MRI/%s.nii',data_dir,q.SubCCIDc{fileindex(s)});
    
    % Copy MEG and MRI data
    copyfile(q.FileNames.meg{fileindex(s)}, destMEG); 
    copyfile(q.FileNames.mri{fileindex(s)}, destMRI); 
    
    % Convert FIF to SPM
    D = spm_eeg_convert(destMEG);
    
    % Coregistration and Forward model
    F = load(q.FileNames.fid{fileindex(s)});
    
    %remove nose points to avoid inclusion in forward model as headshape points and offsetting scalp surface.
    fids = D.fiducials;
    pnts = fids.pnt;
    nose = find(pnts(:,2)>0 & pnts(:,3)<0);
    pnts(nose,:) = [];
    fids.pnt = pnts;
    D = fiducials(D,fids);
    D.save
    
    S = [];
    S.D = D;
    S.fid_mnicoords = 1;
    S.useheadshape=1;
    S.mri = destMRI;
         
    S.use_rhino = 0;
    S.fid.label.nasion = 'Nasion';
    S.fid.label.lpa = 'LPA';
    S.fid.label.rpa = 'RPA';
    S.forward_meg = 'Single Shell';
    S.fid.mnicoords.lpa = F.fid.mni.mm.lpa; 
    S.fid.mnicoords.rpa = F.fid.mni.mm.rpa;
    S.fid.mnicoords.nasion = F.fid.mni.mm.nas;
    
    D=osl_headmodel(S);
    
    D.save
       
end

%% DN 2018 - Check co-registration - preliminary, may change if turn on Rhino.
% Use MATLAB 2015 for below coreg visualisation as does not render properly in later versions. 

% EEG sensor allignment visualisation
 parfor s = 1:length(subs)
     
             D = spm_eeg_load(sprintf('spmeeg_%s',q.SubCCIDc{fileindex(s)}));
             figure('Position', [100 100 1024 1024]);
             spm_eeg_inv_checkdatareg_DN_front_top(D,1,1);
             print(gcf,[check_coreg_directory_EEG, sprintf('Subject_%s', subs{s}), '_check_coreg_front_top', '.bmp'], '-dbmp'); close(gcf)
 end

 parfor s = 1:length(subs)
     
             D = spm_eeg_load(sprintf('spmeeg_%s',q.SubCCIDc{fileindex(s)}));
             figure('Position', [100 100 1024 1024]);
             spm_eeg_inv_checkdatareg_DN_side_oblique(D,1,1);
             print(gcf,[check_coreg_directory_EEG, sprintf('Subject_%s', subs{s}), '_check_coreg_side_oblique', '.bmp'], '-dbmp'); close(gcf) 
 end
 
% MEG sensor allignment visualisation  
 parfor s = 1:length(subs)
      
            D = spm_eeg_load(sprintf('spmeeg_%s',q.SubCCIDc{fileindex(s)}));
            figure('Position', [100 100 1024 1024]);
            if strcmp(subjects{s}, 'CC710261') % no EEG for this subject so indexing for MEG is (D,1,1)
                    spm_eeg_inv_checkdatareg_DN_front_top(D,1,1);
            else
                    spm_eeg_inv_checkdatareg_DN_front_top(D,1,2);
            end
            print(gcf,[check_coreg_directory_MEG, sprintf('Subject_%s', subs{s}), '_check_coreg_front_top', '.bmp'], '-dbmp'); close(gcf)
     
 end
 
 parfor s = 1:length(subs)
      
            D = spm_eeg_load(sprintf('spmeeg_%s',q.SubCCIDc{fileindex(s)}));
            figure('Position', [100 100 1024 1024]);
            if strcmp(subjects{s}, 'CC710261') % no EEG for this subject so indexing for MEG is (D,1,1)
                    spm_eeg_inv_checkdatareg_DN_side_oblique(D,1,1);
            else
                    spm_eeg_inv_checkdatareg_DN_side_oblique(D,1,2);
            end
            print(gcf,[check_coreg_directory_MEG, sprintf('Subject_%s', subs{s}), '_check_coreg_side_oblique', '.bmp'], '-dbmp'); close(gcf)
     
 end


%% Additional preproc steps & Beamformer

p = parcellation(fullfile('fMRI_parcellation_ds8mm.nii.gz')); 
mni_coords = p.template_coordinates;

parfor s = 1:length(subs)
    
    % Remove bad epochs
    D = spm_eeg_load(sprintf('%s/MEG/spmeeg_%s',data_dir,subs{s}));
    D = osl_detect_artefacts(D,'event_significance',0.1); % or D = osl_detect_artefacts(D,'badchannels',false); to only remove bad epochs 
   
    % Downsample the data
    S = [];
    S.D = D;
    S.fsample_new = 200;
    D = spm_eeg_downsample(S);

    % Filter the data
    S = struct;
    S.D = D;
    S.band = 'bandpass';
    S.freq = [1 45];
    D = spm_eeg_filter(S);
    
    % Run Beamforming
    S = struct;
    S.modalities        = {'MEGMAG', 'MEGPLANAR'}; % 'removed MEGMAG and MEGPLANAR' 
    S.timespan          = [0 Inf];
    S.pca_order         = 50;
    S.type              = 'Scalar';
    S.inverse_method    = 'beamform';
    S.prefix            = '';

    D = osl_inverse_model(D,mni_coords,S);
    
    % Select montage  
    D = D.montage('switch',1);
    
    % Save the data
    D.save;

end

%% Delete unused files  

cd(sprintf('%s/MEG',data_dir))

delete(sprintf('%s/MEG/dspmeeg*',data_dir)) % uncomment if want to save space
delete(sprintf('%s/MEG/spmeeg*',data_dir))  % uncomment if want to save space
delete(sprintf('%s/MEG/*.fif',data_dir))    % uncomment if want to save space

delete(gcp('nocreate'))

%% GLEAN

data = strcat(data_dir,'/MEG/fdspmeeg_',subs,'.mat');

% Shuffle the subjects structure before running GLEAN. Added after discussion
% with Nelson J Trujillo_Barreto. Order of subjects can impact
% analysis given can influence temporal elements (problem for you in that
% subjects as per CC numbers are ordered by age and patients are last (alpabetically follow CC with prefix PP)

rng(5) % use this for reproducibility: to get the same randomised struct everytime 
data = data(randperm(numel(data)));

settings = struct;

% Envelope settings 
settings.envelope.overwrite = 0;
settings.envelope.log       = 0;
settings.envelope.fsample   = 20; % we will get the data downsampled to this frequency
settings.envelope.mask      = fullfile(glean_dir,'MNI152_T1_8mm_brain.nii.gz');

% Subspace settings
settings.subspace.overwrite                         = 0;
settings.subspace.normalisation                     = 'none'; % whether we will normalise the data (mean 0, std 1)
settings.subspace.parcellation.file                 = fullfile(glean_dir,'fMRI_parcellation_ds8mm.nii.gz');
settings.subspace.parcellation.orthogonalisation    = 'symmetric'; % method used for leakage correction
settings.subspace.parcellation.method               = 'spatialBasis'; % method used for the parcellation

% Model settings
settings.model.overwrite   = 0;
settings.model.hmm.nstates = nStates; % no. of states
settings.model.hmm.nreps   = 1; % no. of repetitions (we will keep the best according to the free energy)

%% Run GLEAN
GLEAN = glean.setup(glean_name,data,settings);
glean.run(GLEAN)

save(glean_name,'GLEAN')

% compute temporal properties
settings = struct('plot',0);
GLEAN = glean.temporal_stats(GLEAN,settings);

% show the entire state time courses
glean.plot_timecourse(GLEAN)


%% Interrogating the results - set up groups (4 versus updated 5 groupings based on stable or converted MCI) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Using readtable instead of hacky way above to more easily choose groups and exclude age over 88, poor grey matter etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

group_data = readtable('/home/dn01/CamCAN_F/Imaging_analysis_final_related_docs_July_2018/MASTER_ALL_SESSIONS_SPREADSHEET/MASTER_CamCAN_F_data_3_sessions_updated_31st_Dec_2018_test_debug.xlsx', 'TreatAsEmpty', 'NA'); %TreatAsEmpty parameter required to prevent numbers being imported as strings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign groups from master spreadhseet loaded above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CCID = {}; Age = []; Sex = []; MEGID = {}; ACER = []; Group_classic = []; Group_bio_MCI_AD_versus_MCI_stable = [];


Group_classic_temp = group_data{1:end, {'Group_Control_binary_on_initial_screening','Group_Frail_binary_on_initial_screening', 'Group_MCI_binary_on_initial_screening', 'Group_AD_binary_on_initial_screening'}};
for s=1:size(group_data.CC_ID,1)
    if ~strcmp(group_data{s,1},'NA') && group_data{s,229} == 1 % && ~strcmp(group_data{s,1}, 'CC710758') % can add this to skip missing data or exclude specific bad particpants etc. e.g. 'CC710758' omitted from HMM and thus not processsed in beamformer pipeline
        Group_classic(end+1) = find(Group_classic_temp(s,[1:(size(Group_classic_temp,2))]));
        CCID{end+1} = group_data{s,{'CC_ID'}}; % RemeMber need to add {1} after a cell trying to index to extract values as char and not cells (see TF script around megname)
        MEGID{end+1} = group_data{s,{'MEGID'}};
        Age(end+1) = group_data{s, {'AgeOnFrailSession1'}};
        Sex(end+1) = group_data{s, {'Gender_binary'}};
        ACER(end+1) = group_data{s, {'FrailACE_RTotal__100__QCAndCorrectedForWORLDIfAvailable_require'}};
    end
end  

% Oct2018 biomarker and converted updated, now have master columns
% numbered 1-n (e.g. 1 = control, 2 = frail etc), can use below to also
% define classic group above instead of using bianry separate columns as
% using above

for s=1:size(group_data.CC_ID,1)
    if ~strcmp(group_data{s,1},'NA') && group_data{s,229} == 1 % && ~strcmp(group_data{s,1}, 'CC710758') % can add this to skip missing data or exclude specific bad particpants etc. e.g. 'CC710758' omitted from HMM and thus not processsed in beamformer pipeline
        Group_bio_MCI_AD_versus_MCI_stable(end+1) = group_data{s,{'MCI_positives_grouped_with_converters_vs_MCI_negs_and_stable_ve'}};
        CCID{end+1} = group_data{s,{'CC_ID'}}; % Remeber need to add {1} after a cell trying to index to extract values as char and not cells (see TF script around megname)
        MEGID{end+1} = group_data{s,{'MEGID'}};
        Age(end+1) = group_data{s, {'AgeOnFrailSession1'}};
        Sex(end+1) = group_data{s, {'Gender_binary'}};
    end
end 

%%
groupings_classic = Group_classic; % change this to groupings of interest, e.g. biomarker positive patients versus controls etc
groupings_bio_MCI_AD_and_MCI_stable = Group_bio_MCI_AD_versus_MCI_stable;

% Careful, hacky code warning! Need to change below if change groupings above!
Controls = find(groupings_bio_MCI_AD_and_MCI_stable ==1); % Ensure numbering of groups is consistent, e.g. controls = 1, frail = 2 etc
Frail = find(groupings_bio_MCI_AD_and_MCI_stable ==2); % Ensure numbering of groups is consistent, e.g. controls = 1, frail = 2 etc
MCI_stable = find(groupings_bio_MCI_AD_and_MCI_stable ==3); % Ensure numbering of groups is consistent, e.g. controls = 1, frail = 2 etc
MCI_bio_pos_or_convert = find(groupings_bio_MCI_AD_and_MCI_stable ==4); % Ensure numbering of groups is consistent, e.g. controls = 1, frail = 2 etc
AD = find(groupings_bio_MCI_AD_and_MCI_stable ==5); % Ensure numbering of groups is consistent, e.g. controls = 1, frail = 2 etc

groups = {Controls,Frail,MCI_stable,MCI_bio_pos_or_convert,AD}; %will have to change if set up different groups

grpName = {'Controls','Frail','MCI_stable','MCI_bio_pos_or_convert','AD'}; %will have to change if set up differetn groups

% Temporal properties of the states:
% - nOccurrences: the number of occurrences, or state visits.
% - FractionalOccupancy: how much time is state is visited in proportion.
% - MeanLifeTime:: how much time, on average, the state visits last. 
% - MeanIntervalLength: how much time, on average, passes between two consecutive ocurrences of a state.
% - Entropy: the entropy of the state time course (unpredictabillity of flactuations)

% Define measures
mes     = {'nOccurrences','FractionalOccupancy','MeanLifeTime','MeanIntervalLength','Entropy'};
mesName = {'Number of Occurrences','Fractional Occupancy (%)','Mean Life Time (ms)','Mean Interval Length (ms)','Entropy'};
mesCol  = {[0.6 0.8 1],[0 0.4 0.8],[.42 .37 .73],[1 0.6 0.6],[.2 1 0.6]};
cnt = 1;
    
figure
% Loop to get results for all measures - Controls, Frail, MCI and AD
% Loop now sets all figures axis to same value for easier comparison using
% code above. y lim line to get max SD values for each temporal measure

for mesI = 1:length(mes)
    for grp = 1:length(grpName)
    
        clear g
        clear meanG
        clear errG

        eval(sprintf('g = GLEAN.results.temporal_stats.%s.stats(%s,:);',mes{mesI},grpName{grp})); % get results for the selected measure
        g                       = [g]; % add column for age
        g(any(isnan(g),2),:)    = []; % remove NaN values

        % calculate mean and sem
        meanG                   = mean(g(:,1:nStates));
        errG                    = std(g(:,1:nStates))/sqrt(length(g));

        % plot the results
        subplot(length(mes),length(grpName),cnt) % removed subplot(5,4.cnt) so more dynamic if change group number or temporal measures number
        barwitherr(errG,meanG,0.9,'FaceColor',mesCol{mesI})
        ylabel(mesName{mesI},'FontSize',7)
        xlabel('States','FontSize',7)
        set(gca,'FontSize',7)
        xlim([0.5 nStates+0.5])

        strM = sprintf('max(mean(GLEAN.results.temporal_stats.%s.stats));',mes{mesI});
        strS = sprintf('max(std(GLEAN.results.temporal_stats.%s.stats));',mes{mesI}); 
        ylim([0 eval(strM)+(eval(strS))]); %ylim([0 eval(strM)+(eval(strS)*5)]);

        % prepare data for R
        data        = g;
        outputFile  = sprintf('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/CamCAN_F_provisional_analysis_fids_and_headshape_no_EEG_128_no_template/GLEAN/results_for_use_in_R/%s_%s.csv',mes{mesI}, grpName{grp});
        csvwrite(outputFile, data)  
        cnt=cnt+1;

    end
end


%% Sanity check: make sure that FractionalOccupancy < 50% for all subs

figure
bar(max(GLEAN.results.temporal_stats.FractionalOccupancy.stats,[],2),0.9,'FaceColor','m')
ylabel('Max %Occupancy','FontSize',16)
xlabel('Subjects','FontSize',16)
set(gca,'FontSize',16)
xlim([0.5 135]) 

%% Transition probability matrix
Psubject = glean.state_transitions(GLEAN);
Pgroup = zeros(nStates); 
for j = 1:length(Psubject) 
    Psubject{j}(eye(nStates)==1) = 0; 
    % Normalise probabilities such that they sum up to 1
    for k=1:nStates, Psubject{j}(k,:) = Psubject{j}(k,:) / sum(Psubject{j}(k,:)); end
    Pgroup = Pgroup + Psubject{j} / length(Psubject);
end

figure
imagesc(Pgroup); colorbar; axis square
set(gca,'xtick',1:nStates)
set(gca,'ytick',1:nStates)
ylabel('From state','FontSize',16)
xlabel('To state','FontSize',16)
title('Transition Probability Matrix','FontSize',18)
set(gca,'FontSize',15)


%% View spatial maps
settings            = [];
settings.format     = 'nii';
settings.space      = {'parcel'};
GLEAN = glean.pcorr(GLEAN,settings);

cd(sprintf('%s/GLEAN/envelope_L0_R20_F0-Inf/fMRI_parcellation_ds8mm_Ms_Nn/hmm_K%s/results/pcorr/parcel',data_dir, nStates_str))

gunzip('group_pcorr_0-InfHz.nii.gz')
spm_file_split('group_pcorr_0-InfHz.nii')

addpath(genpath('/imaging/dp01/toolboxes/mni2fs_devel/convert/'))


for i= 1:nStates
    figure
    if i <=9
    nii = mni2fs_load_nii(sprintf('group_pcorr_0-InfHz_0000%d.nii',i));
    cmax = max(abs(nii.img(:)));
    cmin = cmax*0.7;
    mni2fs_composite('mnivol',nii,'clims',[cmin cmax],'climstype', 'abs');   
    title(sprintf('STATE %d',i))
    
    else
    nii = mni2fs_load_nii(sprintf('group_pcorr_0-InfHz_000%d.nii',i));
    cmax = max(abs(nii.img(:)));
    cmin = cmax*0.7;
    mni2fs_composite('mnivol',nii,'clims',[cmin cmax],'climstype', 'abs');   
    title(sprintf('STATE %d',i))
    end
end
