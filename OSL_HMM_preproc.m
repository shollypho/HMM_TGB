%% Preprocessing ready for HMM

%% Setup OSL
%
addpath(genpath('/imaging/hp02/software_n_scripts/OSL/osl/osl-core'))
osl_startup
osl_check_installation
% %open osl_debug_log.txt


%% Define paths

addpath('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/scripts')
addpath('/imaging/hp02/TGB/matlab_scripts');
outdir = '/imaging/hp02/TGB/rest_closed/preprocessing/'; % path where data is located
%glean_name = '/imaging/hp02/TGB/rest_closed/hmm/glean/GLEAN.mat'; % name for this GLEAN analysis
%glean_dir = fullfile(osldir,'example_data','glean_example'); % path for GLEAN stuff (masks etc)
check_coreg_directory_EEG = '/imaging/hp02/TGB/MRI/Check_coreg_fiducials_figures_output_EEG/';
check_coreg_directory_MEG = '/imaging/hp02/TGB/MRI/Check_coreg_fiducials_figures_output_MEG/';

cd(outdir)

nStates = 6; % number of states to be used for GLEAN
nStates_str = num2str(nStates);

%% Get data

megPath = '/imaging/hp02/TGB/rest_closed/MF';%'/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/Maxfiltered_all_MEG_tasks_and_resting_state/';
mriPath = '/imaging/hp02/TGB/MRI/';%'/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/aamod_coreg_extended_1_00001_with_fixed_CC710758_and_7_needing_template/';
fidPath = '/imaging/hp02/TGB/MRI/';%'/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/';

all_data = data_setup; % create the cell array containing all the information about the MF/MRI/FID data

%% Convert MF MEG data to SPM format

% loop through sessions and subjects 1-20, 91:end done
for j = 1:length(all_data) % problems in 21, 22, 29, 30
    
    if exist(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))),'file')
        % We don't need to copy existing files
        fprintf('Already imported %s/%s.mat',all_data{j,1},all_data{j,19});
    else
        % Run the conversion to SPM12 format
        cd(outdir)
        mkdir(all_data{j,1});
        cd(all_data{j,1})
        fif_in = fullfile(all_data{j,17},all_data{j,18});
        
        D = osl_import(fif_in,'outfile',sprintf('%s.mat',all_data{j,19}));
        
        
    end
end

%% Perform coregistration
% Two processing steps are carried out in spm_sss. Firstly, we relabel the
% relevant artefact (ECG/EOG) channels for later use and secondly, we run the
% coregistration. Neither of these steps interact with the MEEG data themselves
% so we keep them as a separate stage.
%

% If true, |use_existing| will prevent rerunning the coregstration on any files which already contain a coregistration
use_existing = true;

% Main loop through subjects and sessions
parfor j = 80:99%length(all_data)
    % Load data in from spm_sss, note that any changes in this loop are saved into the same file.
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))));
    
    % Next we re-label the artefact channels
    D = D.chantype(find(strcmp(D.chanlabels,'EEG062')),'EOG');
    D = D.chanlabels(find(strcmp(D.chanlabels,'EEG062')),'VEOG');
    
    D = D.chantype(find(strcmp(D.chanlabels,'EEG061')),'EOG');
    D = D.chanlabels(find(strcmp(D.chanlabels,'EEG061')),'HEOG');
    
    %%% TGB does not have ECG
    %D = D.chantype(find(strcmp(D.chanlabels,'EEG063')),'ECG');
    %D = D.chanlabels(find(strcmp(D.chanlabels,'EEG063')),'ECG');
    
    D.save();
    
    % Skip coreg if needed
    if use_existing && isfield(D,'inv')
        fprintf('Already coregistered %s\n',D.fname);
        continue
    end
    
    % Coregistration is carried out using a call to osl_headmodel. This
    % function takes a struct detailing how the coregistration should be
    % carried out. Importantly, we specify a D object from spm_sss and a
    % strutrual MRI scan from raw_data.
    % Note that useheadshape is set to false, typically we would run this
    % stage with useheadshape set to true. As the MRI scans included in the
    % download have been defaced to ensure participant anonymity the
    % headshape points on the face and nose can cause rotational errors in
    % the coreg. To avoid this we do not include the headshape points and
    % rely only on the fiducials for the coregistration.
    coreg_settings = struct;
    coreg_settings.D = D.fullfile;
    coreg_settings.mri = sprintf('%s/%s/%s',mriPath,all_data{j,1},all_data{j,6}); % individual MRI scans
    coreg_settings.useheadshape = true;
    if strcmp(all_data{j,1},'P9') && all_data{j,4}==2 % do not use head points
        coreg_settings.useheadshape = false;
    end
    coreg_settings.forward_meg = 'Single Shell';
    coreg_settings.use_rhino = false; %currently getting an error with rhino, wait to hear from david
    coreg_settings.fid.label.nasion='Nasion';
    coreg_settings.fid.label.lpa='LPA';
    coreg_settings.fid.label.rpa='RPA';   
    coreg_settings.fid_mnicoords = true;
    coreg_settings.fid.mnicoords.nasion = [all_data{j,8},all_data{j,9},all_data{j,10}];
    coreg_settings.fid.mnicoords.lpa = [all_data{j,11},all_data{j,12},all_data{j,13}];
    coreg_settings.fid.mnicoords.rpa = [all_data{j,14},all_data{j,15},all_data{j,16}];
    
    D = osl_headmodel(coreg_settings);
    
    % ONLY WORKS WITH RHINO - NOT RUNNING RHINO ATM AS NOT WORKING
    % Next we generate and save a summary image so we can check each
    % coregistration has been carried out sucessfully.
    %         h = report.coreg(D);
    %         report.save_figs(h,outdir,D.fname);
    %         close(h);
    
    
end

%% DN 2018 - Check co-registration - preliminary, may change if turn on Rhino.
% Use MATLAB 2015 for below coreg visualisation as does not render properly in later versions.

% EEG sensor allignment visualisation
for j = 1:length(all_data)
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))));
    figure('Position', [100 100 1024 1024]);
    spm_eeg_inv_checkdatareg_DN_front_top(D,1,1);
    print(gcf,[check_coreg_directory_EEG, sprintf('Subject_%s_s%d_rsc%d', all_data{j,1}, all_data{j,4}, all_data{j,5}), '_check_coreg_front_top', '.bmp'], '-dbmp'); close(gcf)
    close all
    
end
%%
% for j = 1:length(all_data)
%
%     D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))));
%     figure('Position', [100 100 1024 1024]);
%     spm_eeg_inv_checkdatareg_DN_side_oblique(D,1,1);
%     print(gcf,[check_coreg_directory_EEG, sprintf('Subject_%s_s%d_rsc%d', all_data{j,1}, all_data{j,4}, all_data{j,5}), '_check_coreg_side_oblique', '.bmp'], '-dbmp'); close(gcf)
%     close all
% end

% MEG sensor allignment visualisation
for j = 1:length(all_data)
    
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))));
    figure('Position', [100 100 1024 1024]);
    if all_data{j,20}==0 % no EEG for this subject so indexing for MEG is (D,1,1)
        spm_eeg_inv_checkdatareg_DN_front_top(D,1,1);
    else
        spm_eeg_inv_checkdatareg_DN_front_top(D,1,2);
    end
    print(gcf,[check_coreg_directory_MEG, sprintf('Subject_%s_s%d_rsc%d', all_data{j,1}, all_data{j,4}, all_data{j,5}), '_check_coreg_front_top', '.bmp'], '-dbmp'); close(gcf)
    close all
end

% for j = 1:length(all_data)
%
%     D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))));
%     figure('Position', [100 100 1024 1024]);
%     if strcmp(all_data{j,1}, 'CC710261') % no EEG for this subject so indexing for MEG is (D,1,1)
%         spm_eeg_inv_checkdatareg_DN_side_oblique(D,1,1);
%     else
%         spm_eeg_inv_checkdatareg_DN_side_oblique(D,1,2);
%     end
%     print(gcf,[check_coreg_directory_MEG, sprintf('Subject_%s_s%d_rsc%d', all_data{j,1}, all_data{j,4}, all_data{j,5}), '_check_coreg_side_oblique', '.bmp'], '-dbmp'); close(gcf)
%     close all
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start MATLAB pool
openparallel = 1;


if openparallel && numel(gcp('nocreate')) == 0
    %if matlabpool('size')==0
    
    cbupool(20);
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



%s = study( fullfile(config.analysisdir,'spm_sss'),0);

% Set output and working directories
outdir = '/imaging/hp02/TGB/rest_closed/preprocessing';%fullfile( config.analysisdir,'spm_sss_processed' );
%wd = getfullpath('tempdir');
%mkdir(wd)
%if ~exist( outdir )
%    mkdir(outdir);
%end



% Main loop through files within the study object
parfor j = 1:length(all_data)
    
    % Load in MEEG data as an SPM object
    Dpar = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s.mat',all_data{j,1},all_data{j,19}))));
    
    
    % Downsample and copy, or just copy
    %    if D.fsample > 250
    Dpar = spm_eeg_downsample(struct('D',Dpar,'fsample_new',250,'prefix','d')); % Note - downsampling cannot be done in-place using prefix='', it just fails
    %delete(sprintf('%s/%s/%s.mat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    %    else
    %        D = D.copy(getfullpath(fullfile(pwd,wd,D.fname))); % Copy into working directory
    %    end
    
    % Apply a 0.1-100Hz passband filter, plus notch filters
    Dpar = osl_filter(Dpar,[0.1 100],'prefix','f');
    Dpar = osl_filter(Dpar,-[49 51],'prefix','f');
    Dpar = osl_filter(Dpar,-[99 101],'prefix','f');
    
    %tidy up intermediate files:
    delete(sprintf('%s/%s/d%s.mat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/d%s.dat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/fd%s.mat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/fd%s.dat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/ffd%s.mat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    delete(sprintf('%s/%s/ffd%s.dat',outdir,all_data{j,1},all_data{j,19})); % delete converted file once downsampling done
    
end
%% detect artefacts
parfor j = 1:length(all_data)
    all_data{j,1}
    % Load in MEEG data as an SPM object
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/fffd%s.mat',all_data{j,1},all_data{j,19}))));
    
    
    % Apply automatric bad segment detection
    [D] = osl_detect_artefacts(D,'badchannels',false, 'badtimes',true);
    D.save;
    %save(sprintf('%s/%s/bad_times.mat', outdir, all_data{j,1}), 'bad_times');
end

%% ICA artefact detection
parfor j = 1:length(all_data)
    all_data{j,1}
    % Load in MEEG data as an SPM object
    Dpar = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/fffd%s.mat',all_data{j,1},all_data{j,19}))));
    
    
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

%% Manual ICA to identify cardiac signal do ica man again up to c10 for eog
for j = 14%:length(all_data)
    %D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/fffd%s.mat',all_data{j,1},all_data{j,19}))));
    %
    %copy(D, fullfile(outdir,sprintf('%s/icaman_fffd%s.mat',all_data{j,1},all_data{j,19})));
    
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/icaman_fffd%s.mat',all_data{j,1},all_data{j,19}))));
    D = osl_africa(D,'do_ident','manual');
    
    D.save;
end



%% Add in ICA for ECG from Rik's ICA script
% p.interpolatebadforICA = 0; % ZZZ HACK - I am not confident that the interpolation is working properly
%
% addpath(genpath('/imaging/hp02/software_n_scripts/eeglab_current/eeglab13_5_4b')) %Put here the path of a version of EEGLAB that has the fileIO toolbox (currently imaging/local version is too old)
% rmpath('/imaging/hp02/software_n_scripts/eeglab_current/eeglab13_5_4b/functions/octavefunc/signal') % Otherwise this overrides firls, causing a fail later.
% addpath('/imaging/hp02/software_n_scripts/dlmcell');
% % %
% ICA.PCA_MEGdim = 60;                    % Number of PCs for ICA - for ADJUST routine, should be set to nEEG channels (70 in my dataset), but for maxfiltered MEG should not be more than 60
% if isfield(p,'interpolatebadforICA') && p.interpolatebadforICA == 1
%     ICA.PCA_EEGdim = 50;                    % ADJUST suggests PCA dim as the same as number of channels, but I find this unstable when badchannels are interpolated, so would suggest reducing to 50 as below. I have modified ADJUST to accept this.
% else
%     ICA.PCA_EEGdim = 70;                %NB: IF YOU CHANGE THIS YOU ALSO NEED TO CHANGE IT AGAIN BELOW
% end
% ICA.FiltPars_orig = [1 45];              % filter bandpass [1 40] %% [0.05 20];
%
% ICA.TemRelZval = 3;                    % Relative temporal threshold in Z-values
% ICA.SpaRelZval = 2;                     % 3 is too strict? (since high topo correlation anyway)
% ICA.Nperm = 0;                              % turn off bootstrapping (takes far too long and not necessary with artefact topographies)
% ICA.VarThr = 0;                               % variance threshold (could be 100/PCA_dim, but suggest 0);
% ICA.Rseed = 1;                                % to make reproducible (not implemented here)
%
% pathstem = '/imaging/hp02/TGB/matlab_scripts/';
% arttopos = load([pathstem, 'MEGArtifactTemplateTopographies']);
%
% pathstem_back = pathstem; %Fix bug where preprocessed ICs are imported from another file
%
% modalities = {'MEGMAG' 'MEGPLANAR'};
% ref_chans = {'EOG061','EOG062'};
%
% for j = 1:length(all_data)
%     all_data{j,1}
%     cd([outdir, '/', all_data{j,1}])
%     files = dir(fullfile(outdir,sprintf('%s/fd%s.mat',all_data{j,1},all_data{j,19})));
%
%     D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/fd%s.mat',all_data{j,1},all_data{j,19}))));
%     D_back = D; %Fix bug where preprocessed ICs are imported from another file
%     refs = []; % Reference signal for correlating with ICs
%     for m = 1:length(modalities)
%         ICA.mod = modalities(m);
%         ICA.refs.spa = {arttopos.HEOG{m}', arttopos.VEOG{m}', arttopos.ECG{m}'};  % Assumes modalities ordered same way!!!
%         chans{m} = find(strcmp(D.chantype,modalities{m}));
%         ICA.d  = D(chans{m},:);
%         ICA.chanlabels = D.chanlabels(chans{m});
%
%         ICA.FiltPars = [ICA.FiltPars_orig D.fsample];
%         %
%         %Load into EEGlab structure to interpolate bad channels
%         EEG = pop_fileio(files.name,'channels',chans{m});
%         eeglocs = [num2cell(1:length(D.sensors('MEG').label))',num2cell(D.sensors('MEG').chanpos),num2cell(D.sensors('MEG').label)];
%         eeglocs = eeglocs(chans{m},:);
%         dlmcell([modalities{m} 'locations.xyz'],eeglocs)
%         EEG = pop_chanedit(EEG, 'load',{[modalities{m} 'locations.xyz'] 'filetype' 'xyz'});
%         if isfield(p,'interpolatebadforICA') && p.interpolatebadforICA == 1
%             [~,tointerpolate] = ismember(D.badchannels,chans{m});
%             EEG = eeg_interp(EEG,nonzeros(tointerpolate));
%             ICA.d = EEG.data;
%             D(chans{m}(nonzeros(tointerpolate)),:) = ICA.d(nonzeros(tointerpolate),:); %Ensure that bad channels to not contribute to the montage - time consuming step
%         end
%         ICA.PCA_dim = ICA.PCA_MEGdim;
%         ICA.TemAbsPval = .05/ICA.PCA_dim;      %0.05;  % Too liberal if not permuted?
%         ICA.SpaAbsPval = .05/ICA.PCA_dim;        %0.05;  % Too liberal if not permuted?
%         %run Rik's ICA with MEG templating
%         [remove{m},weights{m},TraMat{m},temcor{m},spacor{m},varexpl,ICs] = detect_ICA_artefacts(ICA);
%         save(['MEG_ICs_' num2str(m) '_' files(f).name])
%     end
% end

%% Normalise sensor types
redo_pats = [14, 35, 37, 88, 119, 151];
parfor i = 1:6%parfor j=1:length(all_data)
    j = redo_pats(i);
    all_data{j,1}
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/icaman_fffd%s.mat',all_data{j,1},all_data{j,19}))));
    
    copy(D, fullfile(outdir,sprintf('%s/Nicaman_fffd%s.mat',all_data{j,1},all_data{j,19})));
    
    D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/Nicaman_fffd%s.mat',all_data{j,1},all_data{j,19}))));
    
    
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

%% Filtering for different analyses
filt_bw = [1,45;30,100;0.1, 4; 4,8;8,12;12,30];
filt_nom ={'hmm_1_45', 'hmm_gamm', 'aec_delt', 'aec_thet','aec_alph', 'aec_beta'};

parfor i = 1:6%parfor j=1:length(all_data)
    j = redo_pats(i);
    for f = 1:length(filt_nom)
        
        D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/Nicaman_fffd%s.mat',all_data{j,1},all_data{j,19}))));
        all_data{j,1}
        copy(D, fullfile(outdir,sprintf('%s/%s_Nicaman_fffd%s.mat',all_data{j,1},filt_nom{f},all_data{j,19})));
        
        D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s_Nicaman_fffd%s.mat',all_data{j,1},filt_nom{f},all_data{j,19}))));
        D = osl_filter(D,[filt_bw(f,1), filt_bw(f,2)], 'prefix', '');
        
    end
end
%% Beamformer
% Intialise parcellation and study objects
p = parcellation( 'fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm' );

% from david's script:
%p = parcellation(fullfile('fMRI_parcellation_ds8mm.nii.gz'));
mni_coords = p.template_coordinates;
%cbupool(20)
for j = 143%1:length(all_data)
    %j = redo_pats(i);
    %if all_data{j,22} == 1
        for f = 2%length(filt_nom)
            % Load in MEEG data as an SPM object
            D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('%s/%s_Nicaman_fffd%s.mat',all_data{j,1},filt_nom{f},all_data{j,19}))));
            
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
    %end
end
