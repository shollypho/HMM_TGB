addpath(genpath('/imaging/hp02/software_n_scripts/distributionPlot'));

all_data = data_setup; % create the cell array containing all the information about the MF/MRI/FID data

% Define sample rate
sample_rate = 250;

%% Load in results from envelope data

% Meta data
method = 'envelope';
K = 6;

% Find HMM directory
base = fullfile( '/imaging/hp02/TGB/rest_closed/hmm_gamma/', 'hmm_envelope');
mkdir( fullfile( base, 'figures' ) );

% Load in HMM results
load( fullfile(base, sprintf('envelope_HMM_K%s.mat',num2str(K))) );

% Load in run indices
load( fullfile(base, 'envelope_hmm_data.mat'), 'R','B','runlen' );


% Create basepath for saving results
savebase = fullfile( base,'figures','envelope_HMM_K6');

%% Temporal statistics
%
% Here we compute the global temporal statistics from the Amplitude Envelope
% HMM. These are computed per subject and visualised as violin plots

scan_T = [R(1,2) diff(R(:,2))']; % Indexing individual scan sessions
%%%%%%% change subj_T
%subj_T = sum(reshape(scan_T,4,[])); % Indexing individal subjects

% Compute temporal stats

% Fractional Occupancy is the proportion of time spent in each state
FO = getFractionalOccupancy( Gamma, scan_T, 2);
disp('FO');
% Interval Time is the time between subsequent visits to a state
IT = getStateIntervalTimes( Gamma, scan_T, []);
ITmerged = cellfun(@mean,IT);clear IT
disp('IT');
% Life Times (or Dwell Times) is the duration of visits to a state
LT = getStateLifeTimes( Gamma, scan_T, []);
LTmerged = cellfun(@mean,LT); clear LT
disp('LT');

%% Make summary figures
fontsize = 14;

figure('Color','w');
subplot(1,3,1);
distributionPlot(FO,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
title('Fractional Occupancy');xlabel('State');ylabel('Proportion');grid on;
print([savebase '_temporalstats_FO'],'-depsc')

%figure('Color','w');%
subplot(1,3,2);
distributionPlot(LTmerged ./ sample_rate * 1000,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
title('Life Times');xlabel('State');ylabel('Time (ms)');grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
print([savebase '_temporalstats_LT'],'-depsc')

%figure('Color','w');%
subplot(1,3,3);
distributionPlot(ITmerged ./ sample_rate,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
title('Interval Times');xlabel('State');ylabel('Time (secs)');grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
print([savebase '_temporalstats_IT'],'-depsc')

%% Transistion probability matrix
P_raw = hmm.P;
% We normalize, such that we focus on the between state transitions
for k=1:K
    P_raw(k,k) = 0;
    P_raw(k,:) = P_raw(k,:) / sum(P_raw(k,:));
end
disp('HMM-MAR Probability of transition from state i to state j')
P_raw;
figure; imagesc(P_raw);caxis([0 0.6]);colormap('jet');colorbar
%% GLM ANALYSES here


%% Mean Activation Maps
%
% For the envelope HMM, we take the summary of each state directly from the
% observtion models stored in hmm.state, this provides the information for both
% the mean and functional connectivity results.

% load parcellation
parc = parcellation('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm');

% Activation maps are normalised within each state to allow for simple visualisation of the states topology
% net_mean = zeros(39,size(Gamma,2));
%
% for k = 2%:size(Gamma,2)
%
%     net_mean2(:,k) = zscore( diag(hmm.state(k).Omega.Gam_rate) ./ hmm.state(k).Omega.Gam_shape );
%
% end

% for individual brains for each state
net_mean = zeros(39,1);

for k = 1:K%:size(Gamma,2)
    
    net_mean(:,1) = zscore( diag(hmm.state(k).Omega.Gam_rate) ./ hmm.state(k).Omega.Gam_shape );
    
    
    
    
    % visualise state in OSLEYES
    %parc.osleyes(net_mean);
    
    % Optionally save a nifti of the results, these are used to generate the
    % figures in the paper via HCP Workbench
    parc.savenii( net_mean, [savebase '_meanactivations_',num2str(k)]);
end

%% Node Weight Maps

% As with the mean activation maps, the node weights are normalised to aid visualisation
net_nw = zeros(39,size(Gamma,2));
thresh_mean = zeros(size(Gamma,2),1);
for k = 1:size(Gamma,2)
    
    G = hmm.state(k).Omega.Gam_rate ./ hmm.state(k).Omega.Gam_shape;
    G = G - diag(diag(G));
    nw = sum(G,1)' + sum(G,2);
    net_nw(:,k) = zscore(nw);
    
end

% visualise state in OSLEYES
parc.osleyes(net_nw);

% Optionally save a nifti of the results
parc.savenii( net_nw, [savebase '_nodeweights']);


%% Gplotmatrix & MANOVA

% groups
for i = 1:length(data_incl)
    grp{i} = all_data{data_incl(i),2};
end

%gplotmatrix
figure; gplotmatrix(FO, [], grp', [],'+xo');
figure; gplotmatrix(ITmerged, [], grp', [],'+xo');
figure; gplotmatrix(LTmerged, [], grp', [],'+xo');

% MANOVA
[FOd,FOp,FOstats] = manova1(FO,grp)

%% Get data in format for SPSS analysis
cnt=0;
for i = 1:2:length(data_incl)
    cnt = cnt+1;
    % 1: ID
    spss{cnt,1} = all_data{data_incl(i),1}; % ID
    % 2: Con or Pat
    if strcmp(all_data{data_incl(i),2},'Con')
        spss{cnt,2} = 1;
    else
        spss{cnt,2} = 2;
    end
    
    % 3: Con, PSP or bvFTD
    spss{cnt,3} = all_data{data_incl(i),2};
    
    % Out come measures
    if all_data{data_incl(i),7} == 1 % if drug on first session
        for state=1:K
            % 4: FO drug/Placebo
            spss{cnt,state+3} = FO(i,state); % drug
            spss{cnt,state+11}= FO(i+1,state); % placebo
            % 5: LT drug/Placebo
            spss{cnt,state+19} = LTmerged(i,state); % drug
            spss{cnt,state+27} = LTmerged(i+1,state); % placebo
            % 6: IT drug/Placebo
            spss{cnt,state+35} = ITmerged(i,state); % drug
            spss{cnt,state+43} = ITmerged(i+1,state); % placebo
        end
    else
        for state=1:K
            % 4: FO drug/Placebo
            spss{cnt,state+3} = FO(i+1,state); % drug
            spss{cnt,state+11}= FO(i,state); % placebo
            % 5: LT drug/Placebo
            spss{cnt,state+19} = LTmerged(i+1,state); % drug
            spss{cnt,state+27} = LTmerged(i,state); % placebo
            % 6: IT drug/Placebo
            spss{cnt,state+35} = ITmerged(i+1,state); % drug
            spss{cnt,state+43} = ITmerged(i,state); % placebo
        end
    end
end

%% Split into groups for plotting controls vs pats, drug vs placebo

for st = 1:K
    % counters
    con_dr=0;
    con_pl=0;
    pat_dr=0;
    pat_pl=0;
    for i = 1:length(data_incl)
        
        if strcmp(all_data{data_incl(i),2},'Con') && all_data{data_incl(i),7}==1 % Control and drug group
            con_dr = con_dr+ 1;
            FO_group_split(con_dr,st) = FO(i,st);
            IT_group_split(con_dr,st) = ITmerged(i,st);
            LT_group_split(con_dr,st) = LTmerged(i,st);
        elseif strcmp(all_data{data_incl(i),2},'Con') && all_data{data_incl(i),7}==0 % Control and Placebo group
            con_pl = con_pl+ 1;
            FO_group_split(con_pl,st+K) = FO(i,st);
            IT_group_split(con_pl,st+K) = ITmerged(i,st);
            LT_group_split(con_pl,st+K) = LTmerged(i,st);
        elseif (strcmp(all_data{data_incl(i),2},'PSP') || strcmp(all_data{data_incl(i),2},'bvFTD')) && all_data{data_incl(i),7}==1 % Patient and drug group
            pat_dr = pat_dr+ 1;
            FO_group_split(pat_dr,st+(2*K)) = FO(i,st);
            IT_group_split(pat_dr,st+(2*K)) = ITmerged(i,st);
            LT_group_split(pat_dr,st+(2*K)) = LTmerged(i,st);
        elseif (strcmp(all_data{data_incl(i),2},'PSP') || strcmp(all_data{data_incl(i),2},'bvFTD')) && all_data{data_incl(i),7}==0 % Patient and Placebo group
            pat_pl = pat_pl+ 1;
            FO_group_split(pat_pl,st+(3*K)) = FO(i,st);
            IT_group_split(pat_pl,st+(3*K)) = ITmerged(i,st);
            LT_group_split(pat_pl,st+(3*K)) = LTmerged(i,st);
        end
    end
end

% Make group split figures
%% Make summary figures
fontsize = 14;

figure('Color','w');

for st = 1:K
subplot(3,K,st);
distributionPlot({FO_group_split(1:19,st),FO_group_split(1:19,st+K),FO_group_split(1:19,st+(2*K)),FO_group_split(1:19,st+(3*K))},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
grid on;
if st==1
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
end
%print([savebase '_temporalstats_FO'],'-depsc')

subplot(3,K,st+K);
distributionPlot({LT_group_split(1:19,st)./sample_rate*1000,LT_group_split(1:19,st+K)./sample_rate*1000,LT_group_split(1:19,st+(2*K))./sample_rate*1000,LT_group_split(1:19,st+(3*K))./sample_rate*1000},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Life Times');xlabel('Group');ylabel('Time (ms)');
end
grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
%print([savebase '_temporalstats_LT'],'-depsc')
subplot(3,K,st+(2*K));
distributionPlot({LT_group_split(1:19,st)./sample_rate,LT_group_split(1:19,st+K)./sample_rate,LT_group_split(1:19,st+(2*K))./sample_rate,LT_group_split(1:19,st+(3*K))./sample_rate},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Interval Times');xlabel('Group');ylabel('Time (secs)');
end
grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT'],'-depsc')

end


