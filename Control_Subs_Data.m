%% controls details with MRIs

function Con_Sub = Control_Subs_Data()

%%all controls 1:20, not inlcuding sub 7 in car task though for poor behaviour
%% drug sess = which session was drug 'A' placebo
wd = '/imaging/hp02/TGB/rest_closed/';
maknewdirs = 0;

Con_Sub{1}.Name = {'meg16_0250', 'meg16_0261', 'C1', '25947_mprage_nd.nii' };
Con_Sub{1}.Taldir= '/imaging/lh01/FTD2015/MRIs/C1';
Con_Sub{1}.DrugSess = 1;
Con_Sub{1}.Fids = [-1.2 89.2 -14.4; -79.9 -12.1 -51.8; 71.6 -4.6 -62.3] ; %NAS, LPA, RPA
Con_Sub{1}.CarTask = 1;
Con_Sub{1}.RSC = [1,1,1,1];
Con_Sub{1}.EEG = [1,1];

Con_Sub{2}.Name = {'meg17_0016', 'meg17_0036', 'C2', '25945_mprage_nd.nii'};
Con_Sub{2}.Taldir= '/imaging/lh01/FTD2015/MRIs/C2';
Con_Sub{2}.DrugSess = 2;
Con_Sub{2}.Fids = [1.5 74.3 -13.8; -69.8 -1.2 -47.4; 68.9 -6.4 -55.2 ] ; %NAS, LPA, RPA
Con_Sub{2}.CarTask = 1;
Con_Sub{2}.RSC = [1,1,1,1];
Con_Sub{2}.EEG = [1,1];

Con_Sub{3}.Name = {'meg17_0038', 'meg17_0068', 'C3', 'DATA_0006.nii'};
Con_Sub{3}.Taldir= '/imaging/lh01/FTD2015/MRIs/C3';
Con_Sub{3}.DrugSess = 2;
Con_Sub{3}.Fids = [0.7 73.9 -20.8; -66.8 -18.6 -45.6; 70.8 -18.6 -44.3 ] ; %NAS, LPA, RPA
Con_Sub{3}.CarTask = 1;
Con_Sub{3}.RSC = [1,1,1,1];
Con_Sub{3}.EEG = [1,1];

Con_Sub{4}.Name = {'meg17_0030', 'meg17_0067', 'C4', 'sMEGSTRUCTURALS-0005-00001-000192-01.nii'};
Con_Sub{4}.Taldir= '/imaging/lh01/FTD2015/MRIs/C4';
Con_Sub{4}.DrugSess = 2;
Con_Sub{4}.Fids = [0.2 102.5 -0.1; -72.7 15.5 -37.3 ; 69.2 13.0 -41.1 ] ; %NAS, LPA, RPA
Con_Sub{4}.CarTask = 1;
Con_Sub{4}.RSC = [1,0,1,1];
Con_Sub{4}.EEG = [1,1];

Con_Sub{5}.Name = {'meg17_0020', 'meg17_0044', 'C5', '25944_mprage_nd.nii'};
Con_Sub{5}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/25944/20170704_U-ID38292/mprage_nd/';
Con_Sub{5}.DrugSess = 1;
Con_Sub{5}.Fids = [-1.3 79.1 -2.6; -76.2 -6.9 -50.9; 75.0 -11.9 -47.1 ] ; %NAS, LPA, RPA
Con_Sub{5}.CarTask = 1;
Con_Sub{5}.RSC = [1,1,1,1];
Con_Sub{5}.EEG = [1,1];

Con_Sub{6}.Name = {'meg17_0048', 'meg17_0091', 'C6', '25949_mprage_nd.nii'};
Con_Sub{6}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/25949/20170721_U-ID38299/mprage_nd/';
Con_Sub{6}.DrugSess = 1;
Con_Sub{6}.Fids = [4.3 79.8 -8.4; -72.2 -10.7 -40.9; 73.8 -14.9 -46.5 ] ; %NAS, LPA, RPA
Con_Sub{6}.CarTask = 1;
Con_Sub{6}.RSC = [1,0,1,1];
Con_Sub{6}.EEG = [1,0];

%% not including 7 in car task
Con_Sub{7}.Name = {'meg17_0052', 'meg17_0077', 'C7', '17308_mprage_nd.nii'}; %%not including
Con_Sub{7}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/17308/20170707_U-ID38364/mprage_nd/';
Con_Sub{7}.DrugSess = 1;
Con_Sub{7}.Fids = [2.7 80.4 -6.3; -74.8 -1.5 -36.8; 72.6 -11.6 -41.9 ] ; %NAS, LPA, RPA
Con_Sub{7}.CarTask = 0;
Con_Sub{7}.RSC = [1,1,1,1];
Con_Sub{7}.EEG = [1,1];

Con_Sub{8}.Name = {'meg17_0057', 'meg17_0084', 'C8', '21668_mprage_125.nii'};
Con_Sub{8}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/21668/20140306/mprage_125/';
Con_Sub{8}.DrugSess = 2;
Con_Sub{8}.Fids = [3.9 83.4 -0.9; -76.8 -12.3 -45.0 ; 80.8 -13.6 -43.8 ] ; %NAS, LPA, RPA
Con_Sub{8}.CarTask = 1;
Con_Sub{8}.RSC = [1,1,1,1];
Con_Sub{8}.EEG = [1,0];

Con_Sub{9}.Name = {'meg17_0089', 'meg17_0115', 'C9', 'sMEGSTRUCTURAL-0005-00001-000192-01.nii'};
Con_Sub{9}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/CBU/C9/';
Con_Sub{9}.DrugSess = 1;
Con_Sub{9}.Fids = [1.2 74.2 -13.2; -74.9 -12.8 -46.5 ; 70.3 -16.4 -48.9 ] ; %NAS, LPA, RPA
Con_Sub{9}.CarTask = 1;
Con_Sub{9}.RSC = [1,1,1,1];
Con_Sub{9}.EEG = [1,1];

Con_Sub{10}.Name = {'meg17_0081', 'meg17_0098', 'C10', '19909_mprage_nd.nii'};
Con_Sub{10}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/19909/20170721_U-ID38298/mprage_nd/';
Con_Sub{10}.DrugSess = 2;
Con_Sub{10}.Fids = [ 3.6 82.4 -9.5; -79.5 -7.5 -47.0; 75.7 -25.5 -48.4 ] ; %NAS, LPA, RPA
Con_Sub{10}.CarTask = 1;
Con_Sub{10}.RSC = [1,1,1,1];
Con_Sub{10}.EEG = [1,1];

Con_Sub{11}.Name = {'meg17_0096', 'meg17_0101', 'C11', 'single_subj_T1.nii'};
Con_Sub{11}.Taldir = ''; %%no MRI
Con_Sub{11}.DrugSess = 2;
Con_Sub{11}.Fids = [2.2 83.8 -41.6;  -82.9 -18.9 -53.2; 85.6 -17.1 -53.2] ; %NAS, LPA, RPA
Con_Sub{11}.CarTask = 1;
Con_Sub{11}.RSC = [1,1,1,1];
Con_Sub{11}.EEG = [1,1];

Con_Sub{12}.Name = {'meg17_0124', 'meg17_0119', 'C12', '23151_mprage_sag.nii'};
Con_Sub{12}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/23151/20170728/mprage_sag/';
Con_Sub{12}.DrugSess = 2;
Con_Sub{12}.Fids = [0.2 78.9 0.1; -78.1 -12.5 -33.2; 66.6 -9.9 -49.2 ] ; %NAS, LPA, RPA
Con_Sub{12}.CarTask = 1;
Con_Sub{12}.RSC = [1,1,1,1];
Con_Sub{12}.EEG = [1,1];

Con_Sub{13}.Name = {'meg17_0118', 'meg17_0122', 'C13', '19393_mprage_nd.nii'};
Con_Sub{13}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/19393/20170707_U-ID38294/mprage_nd/'; %%;
Con_Sub{13}.DrugSess = 1;
Con_Sub{13}.Fids = [6.6 80.2 -11.6; -79.8 1.2 -48.4; 72.6 -8.3 -53.5 ] ; %NAS, LPA, RPA
Con_Sub{13}.CarTask = 1;
Con_Sub{13}.RSC = [1,1,1,1];
Con_Sub{13}.EEG = [1,1];

Con_Sub{14}.Name = {'meg17_0113', 'meg17_0156', 'C14', '18066_mprage.nii'};
Con_Sub{14}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00373/18066/20140812/mprage/'; %%
Con_Sub{14}.DrugSess = 1;
Con_Sub{14}.Fids = [ -1.2 80.1 -28.2; -82.2 -11.2 -32.9 ; 81.0 -15.9 -35.3 ] ; %NAS, LPA, RPA
Con_Sub{14}.CarTask = 1;
Con_Sub{14}.RSC = [1,1,1,1];
Con_Sub{14}.EEG = [1,1];

Con_Sub{15}.Name = {'meg17_0123', 'meg17_0114', 'C15', '26083_mprage_sag.nii'};
Con_Sub{15}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/26083/20170818/mprage_sag/'; %
Con_Sub{15}.DrugSess = 1;
Con_Sub{15}.Fids = [1.9 74.2 -4.9; -71.6 0.8 -41.0; 75.3 -10.2 -46.6 ] ; %NAS, LPA, RPA
Con_Sub{15}.CarTask = 1;
Con_Sub{15}.RSC = [1,1,1,1];
Con_Sub{15}.EEG = [1,1];

Con_Sub{16}.Name = {'meg17_0126', 'meg17_0147', 'C16', 'DATA_0006.nii'};
Con_Sub{16}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/26485/20180123_U-ID39649/Series_006_Sag_MPRAGE/'; %
Con_Sub{16}.DrugSess = 1;
Con_Sub{16}.Fids = [3.8 73.3 -4.1; -75.3 -9.3 -24.2; 73.5 -17.9 -35.0 ] ; %NAS, LPA, RPA
Con_Sub{16}.CarTask = 1;
Con_Sub{16}.RSC = [1,1,1,1];
Con_Sub{16}.EEG = [1,1];

Con_Sub{17}.Name = {'meg17_0179', 'meg17_0189', 'C17', 'DATA_0006.nii'};
Con_Sub{17}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/26473/20180119_U-ID39625/Series_006_Sag_MPRAGE/'; %
Con_Sub{17}.DrugSess = 1;
Con_Sub{17}.Fids = [-2.6 74.0 -20.5; -73.0 -23.4 -32.9; 73.4 -24.1 -28.8  ] ; %NAS, LPA, RPACon_Sub{19}.Taldir = '/imaging/lh01/FTD2015/MRIs/C19';
Con_Sub{17}.CarTask = 1;
Con_Sub{17}.RSC = [1,0,1,0];
Con_Sub{17}.EEG = [1,1];

Con_Sub{18}.Name = {'meg17_0198', 'meg17_0210', 'C18', 'DATA_0006.nii'};
Con_Sub{18}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/26448/20180111_U-ID39568/Series_006_Sag_MPRAGE/';
Con_Sub{18}.DrugSess = 1;
Con_Sub{18}.Fids = [2.8 86.0 -2.0; -71.5 6.2 -62.2; 74.4 6.2 -49.9  ] ; %NAS, LPA, RPA
Con_Sub{18}.CarTask = 1;
Con_Sub{18}.RSC = [1,1,1,1];
Con_Sub{18}.EEG = [1,1];

Con_Sub{19}.Name = {'meg17_0208', 'meg17_0195', 'C19', 'sMEGStructurals-0005-00001-000192-01.nii'};
Con_Sub{19}.Taldir = '/imaging/lh01/FTD2015/MRIs/C19';
Con_Sub{19}.DrugSess = 2;
Con_Sub{19}.Fids = [-1.4 90.5 -3.9; -75.2 -11.7 -46.6;  75.1 -11.7 -49.2 ] ; %NAS, LPA, RPA
Con_Sub{19}.CarTask = 1;
Con_Sub{19}.RSC = [1,1,1,1];
Con_Sub{19}.EEG = [1,1];

Con_Sub{20}.Name = {'meg17_215', 'meg17_0231', 'C20', 'CBU_MPRAGE_32chn.nii'};
Con_Sub{20}.Taldir = '/imaging/lh01/FTD2015/MRIs/C20';
Con_Sub{20}.DrugSess = 1;
Con_Sub{20}.Fids = [-2.4 90.1 -5.7; -75.9 -3.2 -39.4;  72.3 -2.0 -44.3] ; %NAS, LPA, RPA
Con_Sub{20}.CarTask = 1;
Con_Sub{20}.RSC = [1,1,1,1];
Con_Sub{20}.EEG = [1,1];

%% CAR Task addin MRI dir and MEG dir and other info
dosubs = [1:6 8:20]; %only these sub numbers included in CAR task analyses
for nn = 1:length(dosubs)
    n = dosubs(nn);
    Con_Sub{n}.MRIdir = sprintf('/imaging/lh01/FTD2015/MRIs/%s/',Con_Sub{n}.Name{3} );
    Con_Sub{n}.megpath{1} = sprintf('%s%s/MF/SPM12analysis/', wd, Con_Sub{n}.Name{1});
    Con_Sub{n}.megpath{2} = sprintf('%s%s/MF/SPM12analysis/', wd, Con_Sub{n}.Name{2});
 
end


%% use this to make dir of MRIs
if maknewdirs
    for n = [8:10 12:18] %subnum
        
        filename = char(Con_Sub{n}.Name(4))
        subname = char(Con_Sub{n}.Name(3))
        taldir = Con_Sub{n}.Taldir
        
        mkdir_text = sprintf('mkdir ''%s''',subname)
        cd_text_tal = sprintf('cd ''%s''', taldir)
        copyfile_text = sprintf('copyfile (''%s'', ''/imaging/lh01/FTD2015/MRIs/%s/'')', filename, subname)
        
        cd '/imaging/lh01/FTD2015/MRIs'
        eval(mkdir_text)
        eval(cd_text_tal)
        eval(copyfile_text)
        
    end
    
end





end %% end function