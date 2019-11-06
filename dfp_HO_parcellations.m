%% Extract eigenvariates from Harvard Oxford parcellations (Ece K, 2018)

% Uses random parcellations from HO atlas 25% threshold, removed
% subcortical areas.

% addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906/external/fieldtrip/'));
% addpath('/imaging/local/software/fieldtrip/fieldtrip-latest/external/ctf/') % needed to convert ctf data
% addpath('/imaging/local/software/fieldtrip/fieldtrip-latest/external/mne/') % needed for the elekta data

clear all; clc

impath='/imaging/ek01/atlases/Harvard_Oxford_25/HO_25_1mm/';
outdir='/imaging/ek01/dfp_pilot_EK/graph_theory/HO_data/';
idx=1:98; % parcel ids

subs = {'1002', '1003', '1004', '1005', ...
    '4001', '4002', '6008','4003', '4004',...
    '2001' , '2004','2003'}; %

% create delta v1 files again

datadir = '/imaging/ek01/dfp_pilot_EK/processed_data/'; % MEG dir
subdir = '/rso/';
frontbit = 'pdeMffftsss_';
lastbit = 'delta.mat';
val = 1;  %
tlength=801;
suffix='0_4Hz';

%%

mkdir(outdir)
% 
% uclsubs={'2001','2003','2004'};
% for s=1:length(uclsubs)
%     cd([datadir uclsubs{s} subdir])
%     theData = [datadir uclsubs{s} subdir frontbit lastbit];
%     load(theData);
%     D.sensors.meg.labelorg=D.sensors.meg.labelold;
%     try
%         D.sensors.meg.balance.G1BR.labelorg=D.sensors.meg.balance.G1BR.labelold;
%         D.sensors.meg.balance.G2BR.labelorg=D.sensors.meg.balance.G2BR.labelold;
%         D.sensors.meg.balance.G3BR.labelorg=D.sensors.meg.balance.G3BR.labelold;
%     end
%     try
%         for i=1:numel(D.other.inv{val}.datareg)
%             if strcmp(D.other.inv{val}.datareg.modality,'MEG')
%                 D.other.inv{val}.datareg.sensors.labelorg=D.other.inv{val}.datareg.sensors.labelold;
%                 D.other.inv{val}.datareg.sensors.balance.G1BR.labelorg = D.other.inv{val}.datareg.sensors.balance.G1BR.labelold;
%                 D.other.inv{val}.datareg.sensors.balance.G2BR.labelorg = D.other.inv{val}.datareg.sensors.balance.G2BR.labelold;
%                 D.other.inv{val}.datareg.sensors.balance.G3BR.labelorg = D.other.inv{val}.datareg.sensors.balance.G3BR.labelold;
%             end
%         end
%     end
%     save([datadir uclsubs{s} subdir frontbit lastbit],'D')
% end


% =========================================================================
% Extract the first eigenvariate
% =========================================================================


load('/imaging/ek01/camcan_f_EK/scripts/MEG_source_vertices_MNI.mat');
vertices = int16(vertices);

files=dir(impath); s=1;
for i=3:length(files)
    fnames{s}=files(i).name;s=s+1;
end; clear files s i

for roi=1:max(idx)
    [M XYZ] = spm_read_vols(spm_vol([impath fnames{roi} ]));
    XYZ_mask = XYZ(:,find(M))';
    member = ismember(vertices,XYZ_mask, 'rows');
    XYZ = double(vertices(member,:));
    midpoint(roi,:)=round(mean(XYZ));
    %disp(['Parcel# ' num2str(roi) ' | Vertices# ' num2str(size(XYZ,1))]);
end

warning off;
% (Had to downsample the CTF files with the batch, it was giving errors)
Vall=cell(1,length(subs));

parfor sub=1:length(subs)

    disp([subs{sub} ' ...']); cd([datadir subs{sub} subdir])
    theData = [datadir subs{sub} subdir frontbit lastbit];
    try
        D = spm_eeg_load(theData);
    catch
        D=load(theData);
        D=D.D;
    end
    n=size(D,3)-numel(D.badtrials);
    V=zeros(length(idx),tlength*n);
    D.val = val;
    D.inv{val}.source.fname = '';
    D.inv{val}.source.type  = 'trials'; %'evoked'
    D.inv{val}.source.rad = 0;
            
    for roi=1:length(idx)
        
        disp(['ROI#' num2str(roi)])
        [M, XYZ] = spm_read_vols(spm_vol([impath  fnames{roi}]));
        k=find(M);
        XYZ_mask = XYZ(:,k)';
        member = ismember(vertices,XYZ_mask, 'rows');
        XYZ = double(vertices(member,:)); %clear P M XYZ_mask member k
             
        D.inv{val}.source.XYZ = XYZ;
        Ds = spm_eeg_inv_extract_ac_nosave(D); %removes bad trials
        y=Ds(:,:,:); 
        
        yy=reshape(y,[size(Ds,1),tlength*n]);
        [u, s, u] = svd(yy*yy');
        s = diag(s);
        u       = u(:,1); %first eigenvariate
        v = yy'*u/sqrt(s(1)); 
        %clear t u s v yy
        
        Vall{sub}(roi,:)=v'; %clear Ds rdata       
        
    end %clear XYZ
    
end

for sub=1:length(subs)
    V=squeeze(Vall{sub}(:,:));
    save([outdir 'dfp_' subs{sub} '_HO98_' suffix '.mat'],'V','-v7.3');
    clear V
end
