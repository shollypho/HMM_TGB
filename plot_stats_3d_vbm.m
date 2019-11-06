% plotting parameters/rvalues onto 3D brain parcellations

% using Darren Price's mni2fs toolbox on Github (DPrice80)
% Also help from A Tomassini (see email headed LBA in source space)

%clear all
%close all
clc

% Replace the following path with the path to the mni2fs toolbox folder
toolboxpath = '/imaging/hp02/software_n_scripts/mni2fs-master';
addpath(genpath(toolboxpath)) % will add all subfolders and dependencies

addpath(genpath('/imaging/hp02/software_n_scripts/cbrewer'));
addpath('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/mat_files/source_reconstruction/plotting');
% Harvard Oxford Parcellations
%HO_list % relevant for Ece methods of ROIs
%HO_path = '/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/mat_files/source_reconstruction/Harvard_Oxford';
%NII = load_nii('/imaging/at07/templates/HarvardOxford-cort-maxprob-thr0-2mm.nii');

%numROI = 39;

%% Correlation zvalues


col_gradient = 100;

% min_corr_zval =  -4;%min(pttest.all_action_zero(3,:));
% max_corr_zval = 4;%max(pttest.all_action_zero(3,:));
% range_corr_zval = max_corr_zval- min_corr_zval;
% range_color = range_corr_zval/col_gradient;
% if act_scram == 1
%     prop_color_idx = round((pttest.all_action_zero(3,:)-min_corr_zval)/range_color)+1;
% elseif spec_free == 1;
%     prop_color_idx = round((pttest.all_spec_free(3,:)-min_corr_zval)/range_color)+1;
% end
%
%
% for i = 1:length(prop_color_idx)
%
%     if prop_color_idx(i)<1
%        prop_color_idx(i) = 1;
%     elseif prop_color_idx(i) >col_gradient
%         prop_color_idx(i) = col_gradient;
%     end
%
% end



%% Plotting results


% HO Atlas has 48 ROIs (96 when split l/r)

% colours for the ROIs, will need to think about how to make them for the
% range of ndt or rvalues
c = flipud(cbrewer('seq','YlOrRd',col_gradient));
%c = flipud(cbrewer('seq','Reds',97)); %for ndt

ROIcol = 1:col_gradient;%bins;

%close all
%% Plot the hemispheres
% set up the figure

%close all
% Hemispheres:
hemlab = {'lh','rh'};


inc = 0;


for hem = 1:2
  figure('Color','w','position',[20 72 1100 800])
    hem
    %inc = inc+1;
    
    %ax(inc) = subplot(6,4,inc);
    % Load and Render the FreeSurfer surface
    S = [];
    S.hem = hemlab{hem}; % choose the hemesphere 'lh' or 'rh'
    S.inflationstep = 3;
    S.plotsurf = 'inflated';
    S.lookupsurf = 'mid';
    S.decimation = true; % Decimate the surface for speed.
    S = mni2fs_brain(S);
    
        %
    % overlay the states
       
    %NII = load_nii(sprintf('/imaging/hp02/TGB/rest_closed/hmm_test_retest/hmm_envelope_plc/figures/envelope_HMM_K6_meanactivations_%d.nii',k));
    NII = load_nii('/imaging/hp02/TGB/MRI/Alicia_Preprocessing/c1stats/Classical/spmT_0001.nii');
    
    % Plot an ROI, and make it semi transparent
    S.mnivol = NII;
    S.climstype = 'pos';
    S.colormap = 'hot';%'parula';
    S.clims = [3.38 6];
    S = mni2fs_overlay(S);
    colorbar; caxis([0 6]);
    
    %mni2fs_rescale_colormap
    if hem == 1
        view([-90 0]) % change camera angle
        mni2fs_lights % Dont forget to turn on the lights!
    else
        view([-270 0])
        mni2fs_lights
    end
    
    %saveas(gcf,sprintf('figures/K6overlay_state%d_hem%d_lat.png', k, hem));
    %close all
    
    %inc = inc+1;
    %ax(inc) = subplot(6,4,inc);
%     figure('Color','w','position',[20 72 800 600])
%     S = mni2fs_overlay(S);
%     if hem == 1
%         view([-270 0]) % change camera angle
%         mni2fs_lights % Dont forget to turn on the lights!
%     else
%         view([-90 0])
%         mni2fs_lights
%     end
    % Optional - lighting can be altered after rendering
    %linkaxes(ax)
    %saveas(gcf,sprintf('figures/K6overlay_state%d_hem%d_med.png', k, hem));
   % close all
end

%% Plot colour bar
% figure('Color','w','position',[20 72 400 60]);
% ncol = 96;
% fg=1./ncol; % geometrical factor
%
%  X=fg.*[0 0 1 1];
% Y=0.1.*[1 0 0 1]+(2*iname-1)*0.1;
% for icol=1:ncol
%     X2=X+fg.*(icol-1);
%     fill(X2,Y,c(icol, :), 'linestyle', 'none')
%      xlim([0, 1])
%     hold all
% end %
