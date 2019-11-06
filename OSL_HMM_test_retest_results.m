addpath(genpath('/imaging/hp02/software_n_scripts/distributionPlot'));
addpath(genpath('/imaging/hp02/TGB/matlab_scripts/'));
all_data = data_setup; % create the cell array containing all the information about the MF/MRI/FID data

% Define sample rate
sample_rate = 250;

%% Load in results from envelope data

% Meta data
method = 'envelope';
K = 6;

% Find HMM directory
base = fullfile( '/imaging/hp02/TGB/rest_closed/hmm_test_retest/', 'hmm_envelope_all');
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

% Split Gamma up into  test retest
clear R
R(1,1) = 1;
R(1,2) = 60000;
for i = 2:size(T,1)
    R(i,1)  = R(i-1,1) + 60000;
    R(i,2)  = R(i-1,2) + 60000;
end
for i  = 1:size(T,1)
    R_test1(i,:) = R(i,:);
end
Gamma1 = [];
Gamma2 = [];
for i = 1:2:size(R,1)
    Gamma1 = cat(1,Gamma1,Gamma(R(i,1):R(i,2),:));
end
for i = 2:2:size(R,1)
    Gamma2 = cat(1,Gamma2,Gamma(R(i,1):R(i,2),:));
end


scan_T = [R_test1(1,2) diff(R_test1(:,2))']; % Indexing individual scan sessions

%%%%%%% change subj_T
%subj_T = sum(reshape(scan_T,4,[])); % Indexing individal subjects

% Compute temporal stats

% test1
% Fractional Occupancy is the proportion of time spent in each state
FO_test1 = getFractionalOccupancy( Gamma1, scan_T, 2);
% Interval Time is the time between subsequent visits to a state
IT = getStateIntervalTimes( Gamma1, scan_T, []);
ITmerged_test1 = cellfun(@mean,IT);clear IT
% Life Times (or Dwell Times) is the duration of visits to a state
LT = getStateLifeTimes( Gamma1, scan_T, []);
LTmerged_test1 = cellfun(@mean,LT); clear LT

%test2
% Fractional Occupancy is the proportion of time spent in each state
FO_test2 = getFractionalOccupancy( Gamma2, scan_T, 2);
% Interval Time is the time between subsequent visits to a state
IT = getStateIntervalTimes( Gamma2, scan_T, []);
ITmerged_test2 = cellfun(@mean,IT);clear IT
% Life Times (or Dwell Times) is the duration of visits to a state
LT = getStateLifeTimes( Gamma2, scan_T, []);
LTmerged_test2 = cellfun(@mean,LT); clear LT

% Gamma similarity
[S,assig, gamma1] = getGammaSimilarity (Gamma1, Gamma2);
%% Make summary figures
fontsize = 14;

% test1
figure('Color','w');
subplot(1,3,1);
distributionPlot(FO_test1,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
title('Fractional Occupancy');xlabel('State');ylabel('Proportion');grid on;
%print([savebase '_temporalstats_FO_test1'],'-depsc')
ylim([0 0.5]);

%figure('Color','w');%
subplot(1,3,2);
distributionPlot(LTmerged_test1 ./ sample_rate * 1000,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
title('Life Times');xlabel('State');ylabel('Time (ms)');grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
%print([savebase '_temporalstats_LT_test1'],'-depsc')
ylim([0 500]);

%figure('Color','w');%
subplot(1,3,3);
distributionPlot(ITmerged_test1 ./ sample_rate,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
title('Interval Times');xlabel('State');ylabel('Time (secs)');grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT_test1'],'-depsc')
ylim([0 5]);

% test2

figure('Color','w');
subplot(1,3,1);
distributionPlot(FO_test2,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
title('Fractional Occupancy');xlabel('State');ylabel('Proportion');grid on;
%print([savebase '_temporalstats_FO_test2'],'-depsc')
ylim([0 0.5]);

%figure('Color','w');%
subplot(1,3,2);
distributionPlot(LTmerged_test2 ./ sample_rate * 1000,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
title('Life Times');xlabel('State');ylabel('Time (ms)');grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
%print([savebase '_temporalstats_LT_test2'],'-depsc')
ylim([0 500]);

%figure('Color','w');%
subplot(1,3,3);
distributionPlot(ITmerged_test2 ./ sample_rate,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
title('Interval Times');xlabel('State');ylabel('Time (secs)');grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT_test2'],'-depsc')
ylim([0 5]);

%% GLM ANALYSES here
% Check normal distribution

for i = 1:K
    [FO1_H_norm(i,1), FO1_P_norm(i,1)] = kstest(FO_test1(:,i));
    [FO2_H_norm(i,1), FO2_P_norm(i,1)] = kstest(FO_test2(:,i));
    
%     [IT1_H_norm(i,1), IT1_P_norm(i,1)] = kstest(ITmerged_test1(:,i));
%     [IT2_H_norm(i,1), IT2_P_norm(i,1)] = kstest(ITmerged_test2(:,i));
    
    [LT1_H_norm(i,1), LT1_P_norm(i,1)] = kstest(LTmerged_test1(:,i));
    [LT2_H_norm(i,1), LT2_P_norm(i,1)] = kstest(LTmerged_test2(:,i));
%     
end
% All are non-normal so need to use non-parametric stats below
% for paired ttest use a wilcoxen signed rank test

% Quick stats - paired Wilcoxen
for i = 1:K
    [P_FO(i,1), H_FO(i,1)] = signrank(FO_test1(:,i), FO_test2(:,i), 'tail', 'both');
%     [P_IT(i,1), H_IT(i,1)] = signrank(ITmerged_test1(:,i), ITmerged_test2(:,i), 'tail', 'both');
    [P_LT(i,1), H_LT(i,1)] = signrank(LTmerged_test1(:,i), LTmerged_test2(:,i), 'tail', 'both');
end

% Quick stats: corr
% figure;
for i = 1:K
    [rho_FO(i,1), pcorr_FO(i,1)] = corr(FO_test1(:,i), FO_test2(:,i), 'type', 'Spearman');
    %subplot(3,3,i); scatter(FO_test1(:,i),FO_test2(:,i));
    
%     [rho_IT(i,1), pcorr_IT(i,1)] = corr(ITmerged_test1(:,i), ITmerged_test2(:,i), 'type', 'Spearman');
     [rho_LT(i,1), pcorr_LT(i,1)] = corr(LTmerged_test1(:,i), LTmerged_test2(:,i), 'type', 'Spearman');
end

% ICC
for i = 1:K
[r_FO(i,1), LB_FO(i,1), UB(i,1)] = ICC([FO_test1(:,i), FO_test2(:,i)], 'C-k');
% [r_IT(i,1)] = ICC([ITmerged_test1(:,i), ITmerged_test2(:,i)], 'C-k');
 [r_LT(i,1)] = ICC([LTmerged_test1(:,i), LTmerged_test2(:,i)], 'C-k');
end

% dot plot
% figure; scatter([1:size(FO_test1,1)],FO_test1(:,1)); hold on; scatter([1:size(FO_test2,1)],FO_test2(:,1),'red');

% sortted dots:
[FOsort1,I] = sort(FO_test1(:,2));
FOsort2 = FO_test2(I,2);
% figure; scatter([1:length(FOsort1)],FOsort1); hold on; scatter([1:length(FOsort2)],FOsort2,'red');


% Correlations of state means
for st = 1:K
    mean_FO_t1 = mean(FO_test1);
    mean_FO_t2 = mean(FO_test2);
end

[r_st,p_st] = corr([mean_FO_t1',mean_FO_t2']);

















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


