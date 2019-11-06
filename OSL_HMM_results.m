addpath(genpath('/imaging/hp02/software_n_scripts/'));
addpath(genpath('/imaging/hp02/TGB/matlab_scripts'));
all_data = data_setup; % create the cell array containing all the information about the MF/MRI/FID data

% Define sample rate
sample_rate = 250;

%% Load in results from envelope data

% Meta data
K = 6;

% Find HMM directory
base = fullfile( '/imaging/hp02/TGB/rest_closed/hmm/', 'hmm_envelope', 'run4');
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
subplot(1,2,1);
distributionPlot(FO,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
title('Fractional Occupancy');xlabel('State');ylabel('Proportion');grid on;
%print([savebase '_temporalstats_FO'],'-depsc')

%figure('Color','w');%
subplot(1,2,2);
distributionPlot(LTmerged ./ sample_rate * 1000,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
title('Life Times');xlabel('State');ylabel('Time (ms)');grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
%print([savebase '_temporalstats_LT'],'-depsc')

%figure('Color','w');%
% subplot(1,3,3);
% distributionPlot(ITmerged ./ sample_rate,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
% title('Interval Times');xlabel('State');ylabel('Time (secs)');grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT'],'-depsc')

%% Plot FO and LT with new state order
new_ord = [2 6 5 4 1 3];
FO_neword = zeros(size(FO));
LTmerged_neword = zeros(size(LTmerged));
for i = 1:6
    FO_neword(:,i) = FO(:,new_ord(i));
    LTmerged_neword(:,i) = LTmerged(:,new_ord(i));
end
fontsize = 14;

figure('Color','w');
subplot(1,2,1);
distributionPlot(FO_neword,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
title('Fractional Occupancy');xlabel('State');ylabel('Proportion');grid on;
%print([savebase '_temporalstats_FO'],'-depsc')

%figure('Color','w');%
subplot(1,2,2);
distributionPlot(LTmerged_neword ./ sample_rate * 1000,'showMM',2);%,'color',{set1_cols{1:size(FO,2)}})
title('Life Times');xlabel('State');ylabel('Time (ms)');grid on;

%% Transistion probability matrix
P_raw = hmm.P;
% We normalize, such that we focus on the between state transitions
for k=1:K
    P_raw(k,k) = 0;
    P_raw(k,:) = P_raw(k,:) / sum(P_raw(k,:));
end
disp('HMM-MAR Probability of transition from state i to state j')
P_raw;
figure; imagesc(P_raw);caxis([0 0.5]);colormap('jet');colorbar

P_0thresh=P_raw;
P_0thresh(P_0thresh<0.1) = 0;

% New state order = [2 6 5 4 1 6] or [5 1 6 4 3 2]
P_neword = zeros(6,6);
new_ord = [2 6 5 4 1 3];
for i=1:6
    for j=1:6
        P_neword(i,j) = P_0thresh(new_ord(i), new_ord(j));
    end
end


% HCA - Hierarchical Cluster Analysis
Y = pdist(P_raw);
Z = linkage(Y);
figure; dendrogram(Z);
ordr = [3, 4, 1, 2, 5, 6];

for i = 1:6
    for j=1:6
        P_reordr(i,j) = P_raw(ordr(i),ordr(j));
    end
end
figure; imagesc(P_reordr);caxis([0 0.5]);colormap('jet');colorbar
%xticklabels({'3', '4', '1', '2', '5', '6'});yticklabels({'3', '4', '1', '2', '5', '6'});

% plot as a circle
% addpath(genpath('/imaging/hp02/software_n_scripts/circularGraph-3a7926b'));
% figure;
% myColorMap = lines(length(P_raw));
% myLabel = {'1', '2', '3', '4', '5', '6'};
% P_raw(P_raw <= (1/length(P_raw))) = 0;
% circularGraph(P_raw,'Colormap',myColorMap,'Label',myLabel);
%% GLM ANALYSES here


%% Mean Activation Maps
%
% For the envelope HMM, we take the summary of each state directly from the
% observtion models stored in hmm.state, this provides the information for both
% the mean and functional connectivity results.

% load parcellation
parc = parcellation('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm.nii.gz');

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
[FOd,FOp,FOstats] = manova1(FO,grp);

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

%% Make group split figures - Violin plots

fontsize = 14;

figure('Color','w');
for st = 1:K
subplot(3,K,st);
distributionPlot({FO_group_split(1:18,st),FO_group_split(1:18,st+K),FO_group_split(1:19,st+(2*K)),FO_group_split(1:19,st+(3*K))},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
grid on;
if st==1
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
end
%print([savebase '_temporalstats_FO'],'-depsc')

subplot(3,K,st+K);
distributionPlot({LT_group_split(1:18,st)./sample_rate*1000,LT_group_split(1:18,st+K)./sample_rate*1000,LT_group_split(1:19,st+(2*K))./sample_rate*1000,LT_group_split(1:19,st+(3*K))./sample_rate*1000},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Life Times');xlabel('Group');ylabel('Time (ms)');
end
grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
%print([savebase '_temporalstats_LT'],'-depsc')
subplot(3,K,st+K*2);
distributionPlot({IT_group_split(1:18,st)./sample_rate,IT_group_split(1:18,st+K)./sample_rate,IT_group_split(1:19,st+(2*K))./sample_rate,IT_group_split(1:19,st+(3*K))./sample_rate},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Interval Times');xlabel('Group');ylabel('Time (secs)');
end
grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT'],'-depsc')

end

%% Make group split figures - Bar graphs
figure('Color','w');

for st = 2:4:6
    subplot(2,K,st);
    meanFO = mean([FO_group_split(:,st),FO_group_split(:,st+K),FO_group_split(:,st+(2*K)),FO_group_split(:,st+(3*K))]);
    stderrorF0 = std([FO_group_split(:,st),FO_group_split(:,st+K),FO_group_split(:,st+(2*K)),FO_group_split(:,st+(3*K))])/sqrt(length(FO_group_split(:,st)));
    
    b=bar(1:4,fliplr(meanFO));
    for c = 1:4
        b.FaceColor = 'flat';
        a=(mod(c,2)+1)*3;
        b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
    end
    hold on
    errorbar(1:4,fliplr(meanFO), fliplr(stderrorF0), '.');
    set(gca, 'XTickLabel', ['Pat_P'; 'Pat_D'; 'Con_P'; 'Con_D']);
    
    grid on;
    if st==2
        ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
        ylim([0 0.32])
    else
        ylim([0 0.1])
    end
    
    
    subplot(2,K,st+K);
    meanLT = mean([LT_group_split(:,st)./sample_rate*1000,LT_group_split(:,st+K)./sample_rate*1000,LT_group_split(:,st+(2*K))./sample_rate*1000,LT_group_split(:,st+(3*K))./sample_rate*1000]);%},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
    stderrorLT=std([LT_group_split(:,st)./sample_rate*1000,LT_group_split(:,st+K)./sample_rate*1000,LT_group_split(:,st+(2*K))./sample_rate*1000,LT_group_split(:,st+(3*K))./sample_rate*1000])/sqrt(length(LT_group_split(:,st+(3*K))./sample_rate*1000));
    b=bar(1:4,fliplr(meanLT));
    for c = 1:4
        b.FaceColor = 'flat';
        a=(mod(c,2)+1)*3;
        b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
    end
    hold on
    errorbar(1:4,fliplr(meanLT), fliplr(stderrorLT), '.');
    set(gca, 'XTickLabel', ['Pat_P'; 'Pat_D'; 'Con_P'; 'Con_D']);
    
    
    if st ==2
        title('Life Times');xlabel('Group');ylabel('Time (ms)');
     ylim([0 320])
    else
        ylim([0 250])
    end
    grid on;
    
    
    % subplot(3,K,st+K*2);
    % meanIT = mean([IT_group_split(:,st)./sample_rate,IT_group_split(:,st+K)./sample_rate,IT_group_split(:,st+(2*K))./sample_rate,IT_group_split(:,st+(3*K))./sample_rate]);%},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
    % stderrorIT=std([IT_group_split(:,st)./sample_rate,IT_group_split(:,st+K)./sample_rate,IT_group_split(:,st+(2*K))./sample_rate,IT_group_split(:,st+(3*K))./sample_rate])/sqrt(length(IT_group_split(:,st+(3*K))./sample_rate));
    % b=bar(1:4,fliplr(meanIT));
    % for c = 1:4
    %     b.FaceColor = 'flat';
    %     a=(mod(c,2)+1)*3;
    %     b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
    % end
    % hold on
    % errorbar(1:4,fliplr(meanIT), fliplr(stderrorIT), '.');
    % set(gca, 'XTickLabel', ['Pat_P'; 'Pat_D'; 'Con_P'; 'Con_D']);
    %
    % if st ==1
    %     title('Interval Times');xlabel('Group');ylabel('Time (secs)');
    % end
    % grid on
    
end





%% Split groups up even further, into PSP vs ConPSP, bv vs Conbv & psp vs bv

for st = 1:K
    % counters
    con_dr_psp=0;
    con_pl_psp=0;
    
    con_dr_bv=0;
    con_pl_bv=0;
    
    psp_dr=0;
    psp_pl=0;
    
    bv_dr=0;
    bv_pl=0;
    
    for i = 1:length(data_incl)
        % PSP split
        if all_data{data_incl(i),23} == 2 && all_data{data_incl(i),7}==1 % Control and drug group
            con_dr_psp = con_dr_psp+ 1;
            FO_psp_split(con_dr_psp,st) = FO(i,st);
            IT_psp_split(con_dr_psp,st) = ITmerged(i,st);
            LT_psp_split(con_dr_psp,st) = LTmerged(i,st);
            
        elseif all_data{data_incl(i),23} == 2 && all_data{data_incl(i),7}==0 % Control and Placebo group
            con_pl_psp = con_pl_psp+ 1;
            FO_psp_split(con_pl_psp,st+K) = FO(i,st);
            IT_psp_split(con_pl_psp,st+K) = ITmerged(i,st);
            LT_psp_split(con_pl_psp,st+K) = LTmerged(i,st);
        elseif strcmp(all_data{data_incl(i),2},'PSP') && all_data{data_incl(i),7}==1 % Patient and drug group
            psp_dr = psp_dr+ 1;
            FO_psp_split(psp_dr,st+(2*K)) = FO(i,st);
            IT_psp_split(psp_dr,st+(2*K)) = ITmerged(i,st);
            LT_psp_split(psp_dr,st+(2*K)) = LTmerged(i,st);
        elseif strcmp(all_data{data_incl(i),2},'PSP')  && all_data{data_incl(i),7}==0 % Patient and Placebo group
            psp_pl = psp_pl+ 1;
            FO_psp_split(psp_pl,st+(3*K)) = FO(i,st);
            IT_psp_split(psp_pl,st+(3*K)) = ITmerged(i,st);
            LT_psp_split(psp_pl,st+(3*K)) = LTmerged(i,st);
        end
        
        % bvftd split
        if all_data{data_incl(i),23} == 4 && all_data{data_incl(i),7}==1 % Control and drug group
            con_dr_bv = con_dr_bv+ 1;
            FO_bv_split(con_dr_bv,st) = FO(i,st);
            IT_bv_split(con_dr_bv,st) = ITmerged(i,st);
            LT_bv_split(con_dr_bv,st) = LTmerged(i,st);
        elseif all_data{data_incl(i),23} == 4 && all_data{data_incl(i),7}==0 % Control and Placebo group
            con_pl_bv = con_pl_bv+ 1;
            FO_bv_split(con_pl_bv,st+K) = FO(i,st);
            IT_bv_split(con_pl_bv,st+K) = ITmerged(i,st);
            LT_bv_split(con_pl_bv,st+K) = LTmerged(i,st);
        elseif strcmp(all_data{data_incl(i),2},'bvFTD') && all_data{data_incl(i),7}==1 % Patient and drug group
            bv_dr = bv_dr+ 1;
            FO_bv_split(bv_dr,st+(2*K)) = FO(i,st);
            IT_bv_split(bv_dr,st+(2*K)) = ITmerged(i,st);
            LT_bv_split(bv_dr,st+(2*K)) = LTmerged(i,st);
        elseif  strcmp(all_data{data_incl(i),2},'bvFTD') && all_data{data_incl(i),7}==0 % Patient and Placebo group
            bv_pl = bv_pl+ 1;
            FO_bv_split(bv_pl,st+(3*K)) = FO(i,st);
            IT_bv_split(bv_pl,st+(3*K)) = ITmerged(i,st);
            LT_bv_split(bv_pl,st+(3*K)) = LTmerged(i,st);
        end
        
        
    end
end


%% Make group split figures - Violin plots
addpath('/imaging/hp02/software_n_scripts/mtit');

fontsize = 14;
 
% PSP
figure('Color','w');
for st = 1:K
subplot(3,K,st);
distributionPlot({FO_psp_split(:,st),FO_psp_split(:,st+K),FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K))},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}});
%bar([FO_psp_split(:,st),FO_psp_split(:,st+K),FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K))]');%, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}});

%set(gca,'YLim',[0 1],'FontSize',fontsize)
grid on;
if st==1
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
end
%print([savebase '_temporalstats_FO'],'-depsc')

subplot(3,K,st+K);
distributionPlot({LT_psp_split(:,st)./sample_rate*1000,LT_psp_split(:,st+(K))./sample_rate*1000,LT_psp_split(:,st+(2*K))./sample_rate*1000,LT_psp_split(:,st+(3*K))./sample_rate*1000},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Life Times');xlabel('Group');ylabel('Time (ms)');
end
grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
%print([savebase '_temporalstats_LT'],'-depsc')
subplot(3,K,st+K*2);
distributionPlot({IT_psp_split(:,st)./sample_rate,IT_psp_split(:,st+(K))./sample_rate,IT_psp_split(:,st+(2*K))./sample_rate,IT_psp_split(:,st+(3*K))./sample_rate},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Interval Times');xlabel('Group');ylabel('Time (secs)');
end
grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT'],'-depsc')


end
mtit('PSP vs Control');

% bvFTD
figure('Color','w');
for st = 1:K
subplot(3,K,st);
distributionPlot({FO_bv_split(:,st),FO_bv_split(:,st+K),FO_bv_split(:,st+(2*K)),FO_bv_split(:,st+(3*K))},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
grid on;
if st==1
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
end
%print([savebase '_temporalstats_FO'],'-depsc')

subplot(3,K,st+K);
distributionPlot({LT_bv_split(:,st)./sample_rate*1000,LT_bv_split(:,st+(K))./sample_rate*1000,LT_bv_split(:,st+(2*K))./sample_rate*1000,LT_bv_split(:,st+(3*K))./sample_rate*1000},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Life Times');xlabel('Group');ylabel('Time (ms)');
end
grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
%print([savebase '_temporalstats_LT'],'-depsc')
subplot(3,K,st+K*2);
distributionPlot({IT_bv_split(:,st)./sample_rate,IT_bv_split(:,st+(K))./sample_rate,IT_bv_split(:,st+(2*K))./sample_rate,IT_bv_split(:,st+(3*K))./sample_rate},'showMM',2, 'xNames', {'C_D', 'C_P','P_D','P_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Interval Times');xlabel('Group');ylabel('Time (secs)');
end
grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT'],'-depsc')

end
mtit('bvFTD vs Control');

% psp vs bvFTD
figure('Color','w');
for st = 1:K
subplot(3,K,st);
distributionPlot({FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K)),FO_bv_split(:,st+(2*K)),FO_bv_split(:,st+(3*K))},'showMM',2, 'xNames', {'PSP_D', 'PSP_P','bv_D','bv_P'});%,'color',{set1_cols{1:size(FO,2)}});
%set(gca,'YLim',[0 1],'FontSize',fontsize)
grid on;
if st==1
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
end
%print([savebase '_temporalstats_FO'],'-depsc')

subplot(3,K,st+K);
distributionPlot({LT_psp_split(:,st+(2*K))./sample_rate*1000,LT_psp_split(:,+(3*K))./sample_rate*1000,LT_bv_split(:,st+(2*K))./sample_rate*1000,LT_bv_split(:,st+(3*K))./sample_rate*1000},'showMM',2, 'xNames', {'PSP_D', 'PSP_P','bv_D','bv_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Life Times');xlabel('Group');ylabel('Time (ms)');
end
grid on;
%set(gca,'FontSize',fontsize,'FontSize',fontsize,'YLim',[0 300])
%print([savebase '_temporalstats_LT'],'-depsc')
subplot(3,K,st+K*2);
distributionPlot({IT_psp_split(:,st+(2*K))./sample_rate,IT_psp_split(:,st+(3*K))./sample_rate,IT_bv_split(:,st+(2*K))./sample_rate,IT_bv_split(:,st+(3*K))./sample_rate},'showMM',2, 'xNames', {'PSP_D', 'PSP_P','bv_D','bv_P'});%,'color',{set1_cols{1:size(FO,2)}})
if st ==1
    title('Interval Times');xlabel('Group');ylabel('Time (secs)');
end
grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT'],'-depsc')

end
mtit('PSP vs bvFTD');


%% Make group split figures - Bar graphs:
addpath('/imaging/hp02/software_n_scripts/mtit');

fontsize = 14;

% PSP vs Control
figure('Color','w');
for st = 2:4:K
subplot(2,K,st);
meanFO = mean([FO_psp_split(:,st),FO_psp_split(:,st+K),FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K))]);
stderrorF0 = std([FO_psp_split(:,st),FO_psp_split(:,st+K),FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K))])/sqrt(length([FO_psp_split(:,st),FO_psp_split(:,st+K),FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K))]));
b=bar(1:4,fliplr(meanFO)); 
for c = 1:4
    b.FaceColor = 'flat';
    a=(mod(c,2)+1)*3;
    b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
end
hold on
errorbar(1:4,fliplr(meanFO), fliplr(stderrorF0), '.')
set(gca, 'XTickLabel', ['psp_P'; 'psp_D'; 'Con_P'; 'Con_D'])

%set(gca,'YLim',[0 1],'FontSize',fontsize)
grid on;
if st==2
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
    ylim([0 0.32])
    else
        ylim([0 0.1])
    end
%print([savebase '_temporalstats_FO'],'-depsc')

subplot(2,K,st+K);
meanLT = mean([LT_psp_split(:,st)./sample_rate*1000,LT_psp_split(:,st+(K))./sample_rate*1000,LT_psp_split(:,st+(2*K))./sample_rate*1000,LT_psp_split(:,st+(3*K))./sample_rate*1000]);
stderrorLT = std([LT_psp_split(:,st)./sample_rate*1000,LT_psp_split(:,st+(K))./sample_rate*1000,LT_psp_split(:,st+(2*K))./sample_rate*1000,LT_psp_split(:,st+(3*K))./sample_rate*1000])/sqrt(length([FO_psp_split(:,st),FO_psp_split(:,st+K),FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K))]));
b=bar(1:4,fliplr(meanLT));
for c = 1:4
    b.FaceColor = 'flat';
    a=(mod(c,2)+1)*3;
    b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
end
hold on
errorbar(1:4,fliplr(meanLT), fliplr(stderrorLT), '.')
set(gca, 'XTickLabel', ['psp_P'; 'psp_D'; 'Con_P'; 'Con_D'])


if st ==2
    title('Life Times');xlabel('Group');ylabel('Time (ms)');
    ylim([0 320])
    else
        ylim([0 250])
    end

grid on;

% subplot(3,K,st+K*2);
% meanIT = mean([IT_psp_split(:,st)./sample_rate,IT_psp_split(:,st+(K))./sample_rate,IT_psp_split(:,st+(2*K))./sample_rate,IT_psp_split(:,st+(3*K))./sample_rate]);
% stderrorIT = std([IT_psp_split(:,st)./sample_rate,IT_psp_split(:,st+(K))./sample_rate,IT_psp_split(:,st+(2*K))./sample_rate,IT_psp_split(:,st+(3*K))./sample_rate])/sqrt(length([FO_psp_split(:,st),FO_psp_split(:,st+K),FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K))]));
% b=bar(1:4,fliplr(meanIT));
% for c = 1:4
%     b.FaceColor = 'flat';
%     a=(mod(c,2)+1)*3;
%     b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
% end
% hold on
% errorbar(1:4,fliplr(meanIT), fliplr(stderrorIT), '.')
% set(gca, 'XTickLabel', ['psp_P'; 'psp_D'; 'Con_P'; 'Con_D'])
% 
% 
% 
% if st ==1
%     title('Interval Times');xlabel('Group');ylabel('Time (secs)');
% end
% grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT'],'-depsc')


end
mtit('PSP vs Control');

% bvFTD vs Control
figure('Color','w');
for st = 2:4:K
subplot(2,K,st);
meanFO = mean([FO_bv_split(:,st),FO_bv_split(:,st+K),FO_bv_split(:,st+(2*K)),FO_bv_split(:,st+(3*K))]);
stderrorF0 = std([FO_bv_split(:,st),FO_bv_split(:,st+K),FO_bv_split(:,st+(2*K)),FO_bv_split(:,st+(3*K))])/sqrt(length([FO_bv_split(:,st),FO_bv_split(:,st+K),FO_bv_split(:,st+(2*K)),FO_bv_split(:,st+(3*K))]));
b=bar(1:4,fliplr(meanFO)); 
for c = 1:4
    b.FaceColor = 'flat';
    a=(mod(c,2)+1)*3;
    b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
end
hold on
errorbar(1:4,fliplr(meanFO), fliplr(stderrorF0), '.')
set(gca, 'XTickLabel', ['bvF_P'; 'bvF_D'; 'Con_P'; 'Con_D'])

%set(gca,'YLim',[0 1],'FontSize',fontsize)
grid on;
if st==2
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
    ylim([0 0.32])
    else
        ylim([0 0.1])
    end
%print([savebase '_temporalstats_FO'],'-depsc')

subplot(2,K,st+K);
meanLT = mean([LT_bv_split(:,st)./sample_rate*1000,LT_bv_split(:,st+(K))./sample_rate*1000,LT_bv_split(:,st+(2*K))./sample_rate*1000,LT_bv_split(:,st+(3*K))./sample_rate*1000]);
stderrorLT = std([LT_bv_split(:,st)./sample_rate*1000,LT_bv_split(:,st+(K))./sample_rate*1000,LT_bv_split(:,st+(2*K))./sample_rate*1000,LT_bv_split(:,st+(3*K))./sample_rate*1000])/sqrt(length([FO_bv_split(:,st),FO_bv_split(:,st+K),FO_bv_split(:,st+(2*K)),FO_bv_split(:,st+(3*K))]));
b=bar(1:4,fliplr(meanLT));
for c = 1:4
    b.FaceColor = 'flat';
    a=(mod(c,2)+1)*3;
    b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
end
hold on
errorbar(1:4,fliplr(meanLT), fliplr(stderrorLT), '.')
set(gca, 'XTickLabel', ['bvF_P'; 'bvF_D'; 'Con_P'; 'Con_D'])


if st ==2
    title('Life Times');xlabel('Group');ylabel('Time (ms)');
ylim([0 320])
    else
        ylim([0 250])
    end

grid on;

% subplot(3,K,st+K*2);
% meanIT = mean([IT_bv_split(:,st)./sample_rate,IT_bv_split(:,st+(K))./sample_rate,IT_bv_split(:,st+(2*K))./sample_rate,IT_bv_split(:,st+(3*K))./sample_rate]);
% stderrorIT = std([IT_bv_split(:,st)./sample_rate,IT_bv_split(:,st+(K))./sample_rate,IT_bv_split(:,st+(2*K))./sample_rate,IT_bv_split(:,st+(3*K))./sample_rate])/sqrt(length([FO_bv_split(:,st),FO_bv_split(:,st+K),FO_bv_split(:,st+(2*K)),FO_bv_split(:,st+(3*K))]));
% b=bar(1:4,fliplr(meanIT));
% for c = 1:4
%     b.FaceColor = 'flat';
%     a=(mod(c,2)+1)*3;
%     b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
% end
% hold on
% errorbar(1:4,fliplr(meanIT), fliplr(stderrorIT), '.')
% set(gca, 'XTickLabel', ['bvF_P'; 'bvF_D'; 'Con_P'; 'Con_D'])
% 
% 
% 
% if st ==1
%     title('Interval Times');xlabel('Group');ylabel('Time (secs)');
% end
% grid on
%set(gca,'FontSize',fontsize,'YLim',[0 3]);
%print([savebase '_temporalstats_IT'],'-depsc')


end
mtit('bvFTD vs Control');

% PSP vs bvFTD
figure('Color','w');
for st = 2:4:K
subplot(2,K,st);
meanFO(1:2) = mean([FO_psp_split(:,st+(2*K)),FO_psp_split(:,st+(3*K))]);
meanFO(3:4) = mean([FO_bv_split(:,st+(2*K)),FO_bv_split(:,st+(3*K))]);
stderrorF0(1) = std(FO_psp_split(:,st+(2*K)))/sqrt(length(FO_psp_split(:,st+(2*K))));
stderrorF0(2) = std(FO_psp_split(:,st+(3*K)))/sqrt(length(FO_psp_split(:,st+(3*K))));
stderrorF0(3) = std(FO_bv_split(:,st+(2*K)))/sqrt(length(FO_bv_split(:,st+(2*K))));
stderrorF0(4) = std(FO_bv_split(:,st+(3*K)))/sqrt(length(FO_bv_split(:,st+(3*K))));
b=bar(1:4,fliplr(meanFO)); 
for c = 1:4
    b.FaceColor = 'flat';
    a=(mod(c,2)+1)*3;
    b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
end
hold on
errorbar(1:4,fliplr(meanFO), fliplr(stderrorF0), '.')
set(gca, 'XTickLabel', ['bvF_P'; 'bvF_D'; 'psp_P'; 'psp_D'])

%set(gca,'YLim',[0 1],'FontSize',fontsize)
grid on;
if st==2
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
    ylim([0 0.32])
    else
        ylim([0 0.1])
    end
%print([savebase '_temporalstats_FO'],'-depsc')

subplot(2,K,st+K);
meanLT(1:2) = mean([LT_psp_split(:,st+(2*K))./sample_rate*1000,LT_psp_split(:,st+(3*K))./sample_rate*1000]);
meanLT(3:4) = mean([LT_bv_split(:,st+(2*K))./sample_rate*1000,LT_bv_split(:,st+(3*K))./sample_rate*1000]);
stderrorLT(1) = std(LT_psp_split(:,st+(2*K))./sample_rate*1000)/sqrt(length(LT_psp_split(:,st+(2*K))));
stderrorLT(2) = std(LT_psp_split(:,st+(3*K))./sample_rate*1000)/sqrt(length(LT_psp_split(:,st+(3*K))));
stderrorLT(3) = std(LT_bv_split(:,st+(2*K))./sample_rate*1000)/sqrt(length(LT_bv_split(:,st+(2*K))));
stderrorLT(4) = std(LT_bv_split(:,st+(3*K))./sample_rate*1000)/sqrt(length(LT_bv_split(:,st+(3*K))));


b=bar(1:4,fliplr(meanLT));
for c = 1:4
    b.FaceColor = 'flat';
    a=(mod(c,2)+1)*3;
    b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
end
hold on
errorbar(1:4,fliplr(meanLT), fliplr(stderrorLT), '.')
set(gca, 'XTickLabel', ['bvF_P'; 'bvF_D'; 'psp_P'; 'psp_D'])


if st==2
    ylabel('Proportion');title('Fractional Occupancy');xlabel('Group');
    ylim([0 320])
    else
        ylim([0 250])
    end

grid on;

% subplot(3,K,st+K*2);
% meanIT(1:2) = mean([IT_bv_split(:,st)./sample_rate,IT_bv_split(:,st+(K))./sample_rate]);
% meanIT(3:4) = mean([IT_bv_split(:,st+(2*K))./sample_rate,IT_bv_split(:,st+(3*K))./sample_rate]);
% stderrorIT(1) = std(IT_psp_split(:,st+(2*K))./sample_rate)/sqrt(length(IT_psp_split(:,st+(2*K))));
% stderrorIT(2) = std(IT_psp_split(:,st+(3*K))./sample_rate)/sqrt(length(IT_psp_split(:,st+(3*K))));
% stderrorIT(3) = std(IT_bv_split(:,st+(2*K))./sample_rate)/sqrt(length(IT_bv_split(:,st+(2*K))));
% stderrorIT(4) = std(IT_bv_split(:,st+(3*K))./sample_rate)/sqrt(length(IT_bv_split(:,st+(3*K))));
% 
% b=bar(1:4,fliplr(meanIT));
% for c = 1:4
%     b.FaceColor = 'flat';
%     a=(mod(c,2)+1)*3;
%     b.CData(c,:) = [0.1*a; 0.1*a;0.1*a];
% end
% hold on
% errorbar(1:4,fliplr(meanIT), fliplr(stderrorIT), '.')
% set(gca, 'XTickLabel', ['bvF_P'; 'bvF_D'; 'psp_P'; 'psp_D'])
% 
% 
% 
% if st ==1
%     title('Interval Times');xlabel('Group');ylabel('Time (secs)');
% end
% grid on
% %set(gca,'FontSize',fontsize,'YLim',[0 3]);
% %print([savebase '_temporalstats_IT'],'-depsc')


end
mtit('bvFTD vs PSP');

%% Making group transition circles


t_len = size(data,2)/length(data_incl);
Gamma_round = round(Gamma);



% checks
gamma_sum = sum(Gamma_round,2);
min(gamma_sum);
max(gamma_sum);
sum(gamma_sum);

% Split data into each individual and each session
for i = 1:length(data_incl)
    jstart = (t_len*i)-t_len+1;
    kfin = (t_len*i);
    Gamma_ss{i}.G = Gamma_round(jstart:kfin, 1:6);
    Gamma_ss{i}.G_sum = sum(Gamma_ss{i}.G,2); % Check for when no 
end

% Count transitions

for ss = 1:length(data_incl)
    ss
    Gamma_ss{ss}.transition = zeros(6,6);
    i=1;
    [state_on,j_prev] = max(Gamma_ss{ss}.G(i,:)); % get the first state that is on
    
    for i = 2:size(Gamma_ss{ss}.G,1)
        
        if Gamma_ss{ss}.G_sum(i) ==1 % If there is a clear on state
            
            [state_on,j_now] = max(Gamma_ss{ss}.G(i,:)); % The current state
            
            % Increment the transition from the previous state to the
            % current state:
            Gamma_ss{ss}.transition(j_prev, j_now) = Gamma_ss{ss}.transition(j_prev, j_now)+1;
            
            j_prev = j_now;
            
        end
    end
    
    % Remove diagonal:
    for i = 1:K
        Gamma_ss{ss}.transition(i,i) = 0;
    end
    
    % Turn into percentages along the rows (from)
    sum_row = sum(Gamma_ss{ss}.transition,2);
    Gamma_ss{ss}.transition_prob = zeros(6,6);
    for i = 1:K
        for j = 1:K
            Gamma_ss{ss}.transition_prob(i,j) = Gamma_ss{ss}.transition(i,j)/sum_row(i);
        end
    end
    
end

% Now put into groups:

Gamma_groups.PSP_drug = zeros(6,6);
Gamma_groups.bvFTD_drug = zeros(6,6);
Gamma_groups.Con_drug = zeros(6,6);
Gamma_groups.PSP_placebo = zeros(6,6);
Gamma_groups.bvFTD_placebo = zeros(6,6);
Gamma_groups.Con_placebo = zeros(6,6);


for i = 1:length(data_incl)
    if strcmp(all_data{data_incl(i),2},'Con') && all_data{data_incl(i),7}==1 % Control and drug group
        Gamma_groups.Con_drug = Gamma_groups.Con_drug + Gamma_ss{i}.transition_prob;
    elseif strcmp(all_data{data_incl(i),2},'Con') && all_data{data_incl(i),7}==0 % Control and placebo group
        Gamma_groups.Con_placebo = Gamma_groups.Con_placebo + Gamma_ss{i}.transition_prob;
    elseif strcmp(all_data{data_incl(i),2},'PSP') && all_data{data_incl(i),7}==1 % PSP and drug group
        Gamma_groups.PSP_drug = Gamma_groups.PSP_drug + Gamma_ss{i}.transition_prob;
    elseif strcmp(all_data{data_incl(i),2},'PSP') && all_data{data_incl(i),7}==0 % PSP and placebo group
        Gamma_groups.PSP_placebo = Gamma_groups.PSP_placebo + Gamma_ss{i}.transition_prob;
    elseif strcmp(all_data{data_incl(i),2},'bvFTD') && all_data{data_incl(i),7}==1 % bvFTD and drug group
        Gamma_groups.bvFTD_drug = Gamma_groups.bvFTD_drug + Gamma_ss{i}.transition_prob;
    elseif strcmp(all_data{data_incl(i),2},'bvFTD') && all_data{data_incl(i),7}==0 % bvFTD and placebo group
        Gamma_groups.bvFTD_placebo = Gamma_groups.bvFTD_placebo + Gamma_ss{i}.transition_prob;
    end
end


% Divide by total numbers to get average:
Gamma_groups.PSP_drug = Gamma_groups.PSP_drug/10;
Gamma_groups.bvFTD_drug = Gamma_groups.bvFTD_drug/9;
Gamma_groups.Con_drug = Gamma_groups.Con_drug/19;
Gamma_groups.PSP_placebo = Gamma_groups.PSP_placebo/10;
Gamma_groups.bvFTD_placebo = Gamma_groups.bvFTD_placebo/9;
Gamma_groups.Con_placebo = Gamma_groups.Con_placebo/19;

% Reorder

Gamma_groups.PSP_drug_neword = zeros(6,6);
Gamma_groups.bvFTD_drug_neword = zeros(6,6);
Gamma_groups.Con_drug_neword = zeros(6,6);
Gamma_groups.PSP_placebo_neword = zeros(6,6);
Gamma_groups.bvFTD_placebo_neword = zeros(6,6);
Gamma_groups.Con_placebo_neword = zeros(6,6);

new_ord = [2 6 5 4 1 3];

for i=1:6
    for j=1:6
        Gamma_groups.PSP_drug_neword(i,j) = Gamma_groups.PSP_drug(new_ord(i), new_ord(j));
        Gamma_groups.bvFTD_drug_neword(i,j) = Gamma_groups.bvFTD_drug(new_ord(i), new_ord(j));
        Gamma_groups.Con_drug_neword(i,j) = Gamma_groups.Con_drug(new_ord(i), new_ord(j));
        Gamma_groups.PSP_placebo_neword(i,j) = Gamma_groups.PSP_placebo(new_ord(i), new_ord(j));
        Gamma_groups.bvFTD_placebo_neword(i,j) = Gamma_groups.bvFTD_placebo(new_ord(i), new_ord(j));
        Gamma_groups.Con_placebo_neword(i,j) = Gamma_groups.Con_placebo(new_ord(i), new_ord(j));
    end
end

