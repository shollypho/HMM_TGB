
%% HMM Amplitude envelope
clear all
clc
%% Setup OSL
%
addpath(genpath('/imaging/hp02/software_n_scripts/OSL/osl/osl-core'))
osl_startup
osl_check_installation
% %open osl_debug_log.txt

%% Define paths

addpath('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/scripts')
addpath('/imaging/hp02/TGB/matlab_scripts');
outdir = '/imaging/hp02/TGB/rest_closed/hmm_gamma/hmm_envelope/'; % path where outputs are located

cd(outdir)

nStates = 6; % number of states to be used for hmm

megPath = '/imaging/hp02/TGB/rest_closed/preprocessing/';
mriPath = '/imaging/hp02/TGB/MRI/';

all_data = data_setup; % create the cell array containing all the information about the MF/MRI/FID data

% Find preprocessed data files
%datapath = fullfile( megPath);
%s = study(datapath,'Giles');

%%

% Initialise parcellation:
p = parcellation('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm');

% We'll need to collect the data, T
data = [];            % HMM ready dataset
T = [];               % Length of continuous good segments
R = [];               % Indices of single run within data
B = cell(length(all_data),1);      % Indices of bad samples per session
runlen = zeros(length(all_data),1);          % Length of run per good segment

%%
tic
inc=0;
for j = 1:length(all_data)
    if all_data{j,22} == 1 % inclusion criteria in data_setup function.
        all_data{j,1}
         inc = inc+1;
         data_incl(inc)= j; % keep a record of the indices of the data included
        
        fprintf('Processing %s session %d rsc%d\n',all_data{j,1},all_data{j,4},all_data{j,5});
        
        %-------------------------------
        % continuous file
        D = spm_eeg_load(getfullpath(fullfile(megPath,sprintf('%s/hmm_gamm_Nicaman_fffd%s.mat',all_data{j,1},all_data{j,19}))));
        %D_orig = D.montage('switch',0);
        
        % Do parcellation
        D = ROInets.get_node_tcs(D,p.parcelflag,'spatialBasis','Giles');
    
        
        runlen(j) = size(D,2);
        
        %-------------------------------
        % get power envelope
        dat = osl_envelope( D, 'filter', [30 100], 'orthogonalise',true);
        
        
        
        
        %-------------------------------
        % Get badsamples
        
        
        indx = []; indx = 1:size(dat,2);
        dat = dat(:,good_samples(D));
        indx = indx(:,good_samples(D));
        
              
        
        %------------------------------
        % Remove non-hpi time
        hpi_idx = all_data{j,21};
        % In removing bad samples, may have taken some HPI off time too
        i = 1;
        while indx(i)<=hpi_idx % find the corresponding index
            i = i+1;
        end
        
        
        dat = dat(:,i:end);
        
        %------------------------------
        % Add rsc2 for a couple of participants with short rsc1
        if j == 57 || j == 142 % add envelope for rsc2 for these two participants b/c their rsc1 are too short
            D2 = spm_eeg_load(getfullpath(fullfile(megPath,sprintf('%s/hmm_gamm_Nicaman_fffd%s.mat',all_data{j+1,1},all_data{j+1,19}))));
            D2 = ROInets.get_node_tcs(D2,p.parcelflag,'spatialBasis','Giles');
            dat2 = osl_envelope( D2, 'filter', [2 40], 'orthogonalise',true);
            indx = []; indx = 1:size(dat2,2);
            dat2 = dat2(:,good_samples(D2));
            indx = indx(:,good_samples(D2));
                     
                        %------------------------------
            % Remove non-hpi time
            hpi_idx = all_data{j+1,21};
            % In removing bad samples, may have taken some HPI off time too
            i = 1;
            while indx(i)<=hpi_idx % find the corresponding index
                i = i+1;
            end
            
            
            dat2 = dat2(:,i:end);
            dat = [dat dat2]; % add rsc1 and rsc2 together
        end
       
        datlen(inc,1) =size(dat,2);
        
        %-------------------------------
        % Smooth and normalise
        dat = movmean(dat,25,2,'omitnan'); % 100ms smoothing window
        for ll = 1:p.n_parcels
            dat(ll,:) = ( dat(ll,:) - nanmean(dat(ll,:)) ) ./ nanstd(dat(ll,:));
        end
        %------------------------------
        % Cut everyone to the same length
        dat = dat(:,1:60000);
        t = size(dat,2);
        
        %--------------------------------
        % Store info
        
        offset = sum(T);
        
        R = cat(1,R,[offset+1 offset+size(dat,2)]);
        
        T = cat(1,T,t);
        
        data = cat(2,data,dat);
        
        %--------------------------------
        % Check evoked result
        
        % get trial info
        % %     trl{j} = epochinfo{j}.trl;
        % %
        % %     % replace bad samples with nans
        % %     subj_data = nan(p.n_parcels,size(D,2));
        % %     subj_data(:,good_inds) = data(:,R(j,1):R(j,2));
        % %
        % %     for kk = 1:size(trl{j},1)
        % %         start = trl{j}(kk,1);
        % %         stop = trl{j}(kk,2);
        % %
        % %         if stop > size(subj_data,2) || start > size(subj_data,2)
        % %             disp('some trials missing');
        % %             continue
        % %         else
        % %            erf(:,:,kk,j) = subj_data(:,start:stop)';
        % %         end
        % %     end
   
    end
    
end
toc
% Define HMM folder and save HMM-ready data
hmm_folder = outdir;
if ~exist( hmm_folder )
    mkdir( hmm_folder );
end
outfile = fullfile( hmm_folder, 'envelope_hmm_data' );
%
save( outfile, 'data', 'R', 'T', 'B', 'runlen','data_incl', '-v7.3');

%% HMM inference
%
% Here we infer the HMM itself, a detailed description of the HMM-MAR toolbox
% can be found on https://github.com/OHBA-analysis/HMM-MAR/wiki
%

% Prepare options structure
options = struct();
options.verbose = 1;

% These options specify the data and preprocessing that hmmmar might perform. Further options are discussed here
options.onpower = 0;
options.standardise = 0;
options.Fs = 250;

% Here we specify the HMM parameters
options.K = 6;  	  % The number of states to infer
options.order = 0; 	  % The lag used, this is only relevant when using MAR observations
options.zeromean = 0; 	  % We do want to model the mean, so zeromean is set off
options.covtype = 'full'; % We want to model the full covariance matrix

% These options specify parameters relevant for the Stochastic inference. They
% may be omitted to run a standard inference, but this will greatly increase
% the memory and CPU demands during processing. A detailed description of the
% Stochastic options and their usage can be found here:
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#stochastic

options.BIGNinitbatch = 15;
options.BIGNbatch = 15;
options.BIGtol = 1e-7;
options.BIGcyc = 500;
options.BIGundertol_tostop = 5;
options.BIGdelay = 5;
options.BIGforgetrate = 0.7;
options.BIGbase_weights = 0.9;

% The following loop performs the main HMM inference. We start by
% estimating a 8 state HMM as used in the manuscript.
states_to_infer = [nStates];

% Optionally, we can explore a wider range of values for K by looping through
% several values. This can be done by uncommenting the line below.
% Warning: This is likely to be extremely time-consuming to infer!

%states_to_infer = 2:2:12; % uncomment this line to explore different numbers of states

% The HMM inference is repeated a number of times and the results based on
% the iteration with the lowest free energy. Note that this can be
% extremely time-consuming for large datasets. For a quick exploration of
% results, nrepeats can be set to a smaller value or even 1. The full inference
% is run over 10 repeats.
nrepeats = 1;

% OPEN PARALELL POOL
 cbupool(20)
 
 %% HMM

for kk = states_to_infer
    best_freeenergy = nan;
    options.K = kk;
    
    for irep = 1:nrepeats
        % Run the HMM, note we only store a subset of the outputs
        % more details can be found here: https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#estimation
        tic
        [hmm_iter, Gamma_iter, ~, vpath_iter, ~, ~, ~, ~, fehist] = hmmmar (data',T',options);
        toc
        if isnan(best_freeenergy) || fehist(end) < best_freeenergy
            hmm = hmm_iter;
            Gamma = Gamma_iter;
            vpath = vpath_iter;
        end
    end
    % Save the HMM outputs
    hmm_outfile = fullfile( hmm_folder,  sprintf('envelope_HMM_K%d',options.K));
    save( hmm_outfile ,'hmm','Gamma','vpath','T')
end

