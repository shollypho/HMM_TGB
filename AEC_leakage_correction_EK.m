%% Leakage correction and connectivity script (Ece K, 2018)
%
% This script takes the output from the random parcellations script,
% computes the Hilbert envelopes of each parcel's signal, applies pairwise
% symmetric leakage correction (Colclough et al, 2015), then computes the
% amplitude envelope correlation.

% Code for correction is taken from MEG ROI nets toolbox. The method is
% described as "Apply orthogonalisation on the parcel time-courses. This
% produced orthonormal parcel time-courses which are as close as possible
% to the original time-courses".

% AEC is computed as defined in the Colclough et al 2016 paper: linear
% correlation between the logarithm of ROI power envelopes. We will repeat
% the AEC calculations for the five frequency bands, and form adjacency
% matrices.

%addpath /imaging/ek01/camcam_f_EK/scripts/aal/

clear all; clc;

datadir = '/imaging/ek01/dfp_pilot_EK/graph_theory/HO_data/';
subs={'1002','1003','1004','1005',...
    '4001','4002','4003','4004','6008',...
    '2001','2003','2004'};%
filepart1 = 'dfp_';
filepart2 = '_HO98_';
filepart3 = '0_4Hz';
visit = 'v1';
outdir = '/imaging/ek01/dfp_pilot_EK/graph_theory/HO_matrices/';

%impath='/imaging/ek01/atlases/Craddock_2011/';
idx=1:98; % parcel ids


%%

if ~exist(outdir)
    mkdir(outdir)
end
cd(datadir)

% =====================================================================
% Define parcel locations and new order (run once)
% =====================================================================

% load('/imaging/ek01/camcan_f_EK/scripts/MEG_source_vertices_MNI.mat');
% vertices = int16(vertices); cor_roi=[]; centroid=[];
%
% atlasname = '/imaging/ek01/camcan_f_EK/scripts/aal/ROI_MNI_V5.nii'; % Load AAL atlas to label your clusters
% MNI = spm_read_vols(spm_vol(atlasname));
% [aal1 aal2 aal3] = textread('/imaging/ek01/camcan_f_EK/scripts/aal/ROI_MNI_V5.txt','%s%s%d');
% load('/imaging/ek01/atlases/Craddock_2011/cortical_rois.mat');
%
% clear parcel_label
% for roi=1:80
%     [M XYZ] = spm_read_vols(spm_vol([impath 'Craddock_tcorr05_2level_ROI' num2str(roi) '_r1.nii' ]));
%     XYZ_mask = XYZ(:,find(M))';
%     member = ismember(vertices,XYZ_mask, 'rows');
%     XYZ = mean(double(vertices(member,:)),1);
%     if ~isnan(XYZ(1))
%         mcoor=mni2spmcoor(XYZ);
%         alllabels=nonzeros(MNI(mcoor(:,1)-2:mcoor(:,1)+2,mcoor(:,2)-2:mcoor(:,2)+2,mcoor(:,3)-2:mcoor(:,3)+2)); % search inside a 5mm cube
%         uni_labels=nonzeros(unique(alllabels));
%         if ~isempty(uni_labels)
%             for l=1:length(uni_labels)
%                 perc(l)=length(find(alllabels==uni_labels(l)))/numel(alllabels);
%             end
%             parcel_label{roi}=aal2(find(aal3==uni_labels(find(max(perc)))));
%         else
%             parcel_label{roi}='Outside';
%         end
%     else
%         parcel_label{roi}='No vertices';
%     end
%     clear perc mcoor alllabels uni_labels XYZ
% end
%
% parcel_label=parcel_label';

% manually ordered parcels by lobes

[lobe, area, ~, cor_idx, new_idx] = textread('/imaging/ek01/dfp_pilot_EK/graph_theory/HO_data/parcel_order.txt','%s%s%s%d%d'); % in the order of cor_roi


for sub=1:length(subs)
    disp(['Loading ' subs{sub} ' ...'])
    clear V Vn Vh Vo Vmat
    % =====================================================================
    % Compute Hilbert envelope
    % =====================================================================
    disp(['Computing Hilbert envelope'])
    load([datadir filepart1 subs{sub} filepart2 filepart3 '.mat']); % V
    %V=V(cor_idx,:);
    Vn(new_idx,:)=V(:,:);% put in the order of lobes
    V=Vn;
    
    for p=1:size(V,1)
        Vh(p,:)=abs(hilbert(V(p,:))); % parcel analytic signal
    end
    
    
    % =====================================================================
    % Apply symmetric & pairwise univariate leakage correction and compute AEC
    % =====================================================================
    
    % Find the closest orthonormal matrix,L, orthonormal matrix, W,
    % weighting matrix such that L = V * W
    
    %     for i=1:size(Vh,1)
    %         for j=1:size(Vh,1)
    %             if ~isequal(i,j)
    %                 [Vo, ~, ~, W] = symmetric_orthogonalise(Vh([i,j],:)',true);
    %                 Vo=Vo';
    %                 clean = all(isfinite(Vo),1);
    %                 c=corr(Vo(:,clean).');
    %                 Vmat(i,j) = abs(c(1,2));
    %             else
    %                 Vmat(i,j)=1;
    %             end
    %         end
    %     end
    
    % =====================================================================
    % Apply symmetric & multivariate leakage correction and compute AEC
    % =====================================================================
    
    disp('Correcting leakage and computing AEC')
    [Vo, ~, ~, W] = symmetric_orthogonalise(Vh',true);
    Vo=Vo';
    clean = all(isfinite(Vo),1);
    
    Vmat=corr(log10(abs(Vo(:,clean))).'); %normalised Hipp et al
    %Vmat2=corr(Vo(:,clean).'); %raw
   
    all_Vmat(sub,:,:)=Vmat;
    %all_Vmat2(sub,:,:)=Vmat2;
    
    save([outdir filepart1 subs{sub} filepart2 strrep(filepart3,'_v2','') '_' visit '_asso.mat'], 'Vmat');
    disp('Done!')
end


% =========================================================================
% Plot figures
% =========================================================================

mean_Vmat=squeeze(mean(all_Vmat,1));
%mean_Vmat2=squeeze(mean(all_Vmat2,1));
ticks=[11 29 45 62 83];
split=[1 22; 23 36; 37 54; 55 68; 69 98];
lobes={'FRO','LIM','OCC','PAR','TEM'};

figure; imagesc(mean_Vmat); 
colormap('jet'); colorbar; caxis([0 0.1]);daspect([1 1 1]); 
title(['Connectivity mean in ' strrep(filepart3,'_','-')]); 
set(gca,'XTick',ticks, 'XTickLabel',lobes); 
set(gca,'YTick',ticks, 'YTickLabel',lobes); hold on;

for i=1: size(split,1)
    rectangle('Position',[split(i,1) split(i,1) split(i,2)-split(i,1) split(i,2)-split(i,1)], ...
        'LineWidth',3,'EdgeColor','b')
end

print(gcf,[outdir 'group_cormat_norm_' filepart3 '_' visit '.bmp'], '-dbmp'); close(gcf)

% figure; imagesc(mean_Vmat2); 
% colormap('jet'); colorbar; caxis([-0.05 0.05]); daspect([1 1 1]); title(filepart3)
% set(gca,'XTick',split(1:end-1), 'XTickLabel',lobes); xtickangle(90)
% set(gca,'YTick',split(1:end-1), 'YTickLabel',lobes); hold on;
% for l=2:length(split)-1
%     plot([0.5 size(mean_Vmat,2)+0.5],[split(l)+0.5 split(l)+0.5],'k','LineWidth',2)
%     plot([split(l)+0.5 split(l)+0.5], [0.5 size(mean_Vmat,2)+0.5],'k','LineWidth',2)
% end
% print(gcf,[outdir 'group_cormat_raw_' filepart3 '.bmp'], '-dbmp'); close(gcf)
% hold on; 





