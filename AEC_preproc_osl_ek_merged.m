%% Preprocessing ready for AEC
% Merging OSL pipeline with Ece's AEC pipeline because these might be in
% the same paper to best to keep the preprocessing as similar as possible.


%% Setup OSL
%
 addpath(genpath('/imaging/hp02/software_n_scripts/OSL/osl/osl-core'))
 osl_startup
 osl_check_installation
% %open osl_debug_log.txt


%% Define paths

addpath('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/scripts')
addpath('/imaging/hp02/TGB/matlab_scripts');
outdir = '/imaging/hp02/TGB/rest_closed/aec/'; % path where data is located
%glean_name = '/imaging/hp02/TGB/rest_closed/hmm/glean/GLEAN.mat'; % name for this GLEAN analysis
%glean_dir = fullfile(osldir,'example_data','glean_example'); % path for GLEAN stuff (masks etc)
check_coreg_directory_EEG = '/imaging/hp02/TGB/MRI/Check_coreg_fiducials_figures_output_EEG/';
check_coreg_directory_MEG = '/imaging/hp02/TGB/MRI/Check_coreg_fiducials_figures_output_MEG/';

cd(outdir)

nStates = 8; % number of states to be used for GLEAN
nStates_str = num2str(nStates);

%% Get data

megPath = '/imaging/hp02/TGB/rest_closed/MF';%'/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/Maxfiltered_all_MEG_tasks_and_resting_state/';
mriPath = '/imaging/hp02/TGB/MRI/';%'/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/aamod_coreg_extended_1_00001_with_fixed_CC710758_and_7_needing_template/';
fidPath = '/imaging/hp02/TGB/MRI/';%'/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/';

all_data = data_setup; % create the cell array containing all the information about the MF/MRI/FID data

%% Copy HMM downsampled data, ready for ICA and filtering

% loop through sessions and subjects
parfor j = 1:length(all_data)
    if exist(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))),'file')
        % We don't need to copy existing files
        fprintf('Already imported %s/%s.mat',all_data{j,1},all_data{j,19});
    else
        % Run the conversion to SPM12 format
        cd(outdir)
        mkdir(all_data{j,1});
        cd(all_data{j,1})
        copyfile(sprintf('/imaging/hp02/TGB/rest_closed/hmm/%s/%s.mat',all_data{j,1},all_data{j,19}),'./')
    end
    
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start MATLAB pool
openparallel = 1;


 if openparallel && numel(gcp('nocreate')) == 0
     %if matlabpool('size')==0
         
         cbupool(167);
%          P.ResourceTemplate = '-l nodes=^N^,mem=4GB,walltime=140:00:00';
%          P.SubmitArguments='-W x="NODESET:ONEOF:FEATURES:MAXFILTER"';
%          poolhandle = parpool(P);
     %end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% MEEG Data Preprocessing
%
% All the processing including the electrophysiological data is carried out
% within this loop. Here we apply a range of data preprocessing steps to remove
% artefacts from our data, project the data into source-space and apply a
% parcellation.
%
% Downsampling - Downsample from 1000Hz to 250Hz
% Filtering - We apply a single passband filter to isolate data between 1 and 45Hz
% Bad Segment Detection - An automatic algorithm is used to identify noisy data segments which are removed from subsequent analysis
% Independant Components Analysis - ICA is used to identify artefactual componets by correlation with the EOG and ECG, these are removed from subsequent analysis
% Sensor Normalisation - The Magnetometers and Gradiometers within each dataset are normalised to make their variances comparable.
% Beamforming - An LCMV Beamformer is used to project the sensor data into an 8mm grid in source space
% Parcellation - The source space data is reduced into a network of parcels defined by a NIFTI file
%
% The analysis is carried out on files copied across from spm_sss.  Where
% possible, the processsing is carried out on the files in spm_sss_processed
% 'in place', meaning that they do not generate a new SPM object. The
% downsampling and filtering overwrite the old data file, but the ICA,
% beamforming and parcellation are applied by adding online montages to the SPM
% object.



% Main loop through files within the study object
parfor j = 1:length(all_data)
    
    % Load in MEEG data as an SPM object
    Dpar = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))));
    
    % Downsample and copy, or just copy
    %    if D.fsample > 250
    Dpar = spm_eeg_downsample(struct('D',Dpar,'fsample_new',250,'prefix','d')); % Note - downsampling cannot be done in-place using prefix='', it just fails
    

    % Apply a 0.1-100Hz passband filter
    Dpar = osl_filter(Dpar,[0.1 100],'prefix','f');
    delete(sprintf('%s/%s/d%s.dat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/%s.dat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/d%s.mat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/%s.mat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    
    % Apply notch filter:
    Dpar = osl_filter(Dpar,-[48 52],'prefix', 'f') 
    delete(sprintf('%s/%s/fd%s.mat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/fd%s.dat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    
    all_data{j,1}
    
end
%% detect artefacts
parfor j = 1:length(all_data)
    all_data{j,1}
    % Load in MEEG data as an SPM object
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/ffd%s.mat',all_data{j,1},all_data{j,19}))));
    
    
    % Apply automatric bad segment detection
    [D] = osl_detect_artefacts(D,'badchannels',false, 'badtimes',true);
    D.save;
    %save(sprintf('%s/%s/bad_times.mat', outdir, all_data{j,1}), 'bad_times');
end

%% ICA artefact detection
parfor j = 1:length(all_data)
    all_data{j,1}
    % Load in MEEG data as an SPM object
    Dpar = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/ffd%s.mat',all_data{j,1},all_data{j,19}))));
    
        
    % Run ICA artefact detection. This will automatically reject components
    % which have correlations larger than .5 with either of the artefact
    % channels.
    Dpar = osl_africa(Dpar,'used_maxfilter',true);

end

%% check ICAs

% Preallocate array for ICA sessions that might need manual checking
sessions_to_check = [];

for j = 1:length(all_data)
    % Though the automatic correlations generally work well, we should be
    % careful to check for unsusual datasets and possibly manually correct the
    % automatic assessment. This is particularly important for relatively noisy
    % data or when analysing a new dataset for the first time.
    %
    % Here we will manally inspect the ICA component rejections for any dataset meeting the following criteria
    % # More than 4 ICs rejected
    % # Zero ICs rejected
    % # No component rejected due to correlation with EOG
    % # No component rejected due to correlation with ECG
    
    % Load in MEEG data as an SPM object
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/fd%s.mat',all_data{j,1},all_data{j,19}))));
       
    artefact_chan_corr_thresh = .5;
    if isempty(D.ica.bad_components) || length(D.ica.bad_components) > 4
        sprintf('%s components rejected, recommend checking session', length(D.ica.bad_components))
        sessions_to_check = cat(1,sessions_to_check,j);
    elseif max(D.ica.metrics.corr_chan_367_EOG.value) < artefact_chan_corr_thresh || ...
            max(D.ica.metrics.corr_chan_368_EOG.value) < artefact_chan_corr_thresh
        disp('no candidate components for either EOG, recommend checking session');
        sessions_to_check = cat(1,sessions_to_check,j);
        %elseif max(D.ica.metrics.corr_chan_369_ECG.value) < artefact_chan_corr_thresh
        %          disp('no candidate components for  ECG, recommend checking session');
        %    sessions_to_check = cat(1,sessions_to_check,j);
    end
    
    % This is where manual Africa would go
    % D = D.montage('remove',1:D.montage('getnumber'));
    % D = osl_africa(D,'do_ident','manual');
end
%% This is as far as I've got with AEC for now, need Ece's help on source loc and parcellations (want more than 39 parcels!)


%% Manual ICA to identify cardiac signal
for j = 58:length(all_data)
    if all_data{j,5} == 1 || j==58 || j==143 % save some time and only incl rsc1 and the couple of rsc2s needed to make up for short rsc1s
     D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/ffd%s.mat',all_data{j,1},all_data{j,19}))));
%     
    copy(D, fullfile(outdir,sprintf('%s/icaman_ffd%s.mat',all_data{j,1},all_data{j,19})));
   
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/icaman_ffd%s.mat',all_data{j,1},all_data{j,19}))));
    D = osl_africa(D,'do_ident','manual');
    
    D.save;
    end
end

%% Normalise sensor types
parfor j = 1:length(all_data)
    all_data{j,1}
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/icaman_ffd%s.mat',all_data{j,1},all_data{j,19}))));
    
    copy(D, fullfile(outdir,sprintf('%s/Nicaman_ffd%s.mat',all_data{j,1},all_data{j,19})));
    
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/Nicaman_ffd%s.mat',all_data{j,1},all_data{j,19}))));
    
    
    S = [];
    S.D = D;
    S.modalities = {'MEGMAG','MEGPLANAR'};
    S.do_plots = 0;
    S.samples2use = good_samples(D,D.indchantype(S.modalities,'GOOD'));
    S.trials = 1;
    S.pca_dim = 99;
    S.force_pca_dim = 0;
    S.normalise_method = 'min_eig';
    D = normalise_sensor_data( S );
end

%% Beamformer
% Intialise parcellation and study objects
p = parcellation( 'fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm' );

% from david's script:
%p = parcellation(fullfile('fMRI_parcellation_ds8mm.nii.gz'));
mni_coords = p.template_coordinates;

parfor j = 1:length(all_data)
    % Load in MEEG data as an SPM object
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/Nicaman_fd%s.mat',all_data{j,1},all_data{j,19}))));
        
    % Run LCMV Beamformer
%     S = struct;
%     S.modalities        = {'MEGMAG', 'MEGPLANAR'}; % 'removed MEGMAG and MEGPLANAR'
%     S.timespan          = [0 Inf];
%     S.pca_order         = 50;
%     S.type              = 'Scalar';
%     S.inverse_method    = 'beamform';
%     D = osl_inverse_model(D,p.template_coordinates,S);
    D = osl_inverse_model(D,p.template_coordinates,'pca_order',50);
    
    % Do parcellation
%     D = ROInets.get_node_tcs(D,p.parcelflag,'spatialBasis','Giles');
    
    % Save out
%     D = D.montage('switch',0);
%     D.copy(fullfile(outdir,sprintf('%s/beam_fd%s.mat',all_data{j,1},all_data{j,19})));
    
end
