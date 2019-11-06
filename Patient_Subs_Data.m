%% controls details with MRIs

function Pat_Sub = Patient_Subs_Data()

%%al;l controls 1:20, not inlcuding sub 7 though for poor behaviour
%% drug sess = which session was drug 'A' placebo
wd = '/imaging/lh01/FTD2015/MEG_DATA/ActionTask/';


Pat_Sub{1}.Name = {'NA', 'NA', 'P1', '' };
Pat_Sub{1}.Taldir= '';
Pat_Sub{1}.DrugSess = 0;
Pat_Sub{1}.Fids = [] ; %NAS, LPA, RPA
Pat_Sub{1}.Dia = '';
Pat_Sub{1}.CarTask = 0;
Pat_Sub{1}.RSC= [0,0,0,0];
Pat_Sub{1}.EEG = [];

Pat_Sub{2}.Name = {'meg16_0284', 'meg16_0295', 'P2', '25335_mprage_sag.nii'};
Pat_Sub{2}.Taldir= '/imaging/na01/HPHI_MRIs/tallie/p00415/25335/20161115/mprage_sag/';
Pat_Sub{2}.DrugSess = 1;
Pat_Sub{2}.Fids = [-1.2 80.8 -4.2; -72.3 1.5 -52.2; 69.9 -2.6 -41.2] ; %NAS, LPA, RPA
Pat_Sub{2}.Dia = 'PSP';
Pat_Sub{2}.CarTask = 1;
Pat_Sub{2}.RSC= [1,1,1,1];
Pat_Sub{2}.EEG = [1,1];

Pat_Sub{3}.Name = {'meg16_0287', 'meg16_0308', 'P3', '25584_mprage_sag.nii'};
Pat_Sub{3}.Taldir='/imaging/na01/HPHI_MRIs/tallie/p00415/25584/20161116/mprage_sag/';
Pat_Sub{3}.DrugSess = 2;
Pat_Sub{3}.Fids = [3.6 80.1 -7.9; -73.5 -4.4 -48.9 ; 69.1 -5.7 -63.0 ] ; %NAS, LPA, RPA
Pat_Sub{3}.Dia = 'PSP';
Pat_Sub{3}.CarTask = 1;
Pat_Sub{3}.RSC = [1,0,1,1];
Pat_Sub{3}.EEG = [1,1];

Pat_Sub{4}.Name = {'meg16_0335', 'NA', 'P4', ''};
Pat_Sub{4}.Taldir= '';
Pat_Sub{4}.DrugSess = 2;
Pat_Sub{4}.Fids = [] ; %NAS, LPA, RPA
Pat_Sub{4}.Dia = 'PSP';
Pat_Sub{4}.CarTask = 0;
Pat_Sub{4}.RSC = [0,0,0,0];
Pat_Sub{4}.EEG = [];

Pat_Sub{5}.Name = {'meg17_0002', 'meg17_0060', 'P5', '23936_mprage.nii'};
Pat_Sub{5}.Taldir ='/imaging/na01/HPHI_MRIs/tallie/p00373/23936/20150716/mprage/';                     
Pat_Sub{5}.DrugSess = 2;
Pat_Sub{5}.Fids = [2.8 88.0 5.8; -67.6 0.8 -38.2; 76.8 -9.9 -34.6] ; %HOLLY NAS, LPA, RPA
Pat_Sub{5}.Dia = 'PSP';
Pat_Sub{5}.CarTask = 0;
Pat_Sub{5}.RSC = [1,0,1,1];
Pat_Sub{5}.EEG = [1,1];

Pat_Sub{6}.Name = {'meg17_0005', 'meg17_0028', 'P6', '23503_mprage.nii'};
Pat_Sub{6}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00373/23503/20150626/mprage/';                      
Pat_Sub{6}.DrugSess = 1;
Pat_Sub{6}.Fids = [-1.0 78.2 -15.6; -68.5 -2.9 -37.4; 65.4 -6.4 -50.0 ] ; %NAS, LPA, RPA
Pat_Sub{6}.Dia = 'PSP';
Pat_Sub{6}.CarTask = 1;
Pat_Sub{6}.RSC = [0,0,0,0];% not psp
Pat_Sub{6}.EEG = [];

Pat_Sub{7}.Name = {'meg17_0007', 'meg17_0039', 'P7', '25556_mprage_sag.nii'}; %%not including
Pat_Sub{7}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00259/25556/20161005/mprage_sag/';              
Pat_Sub{7}.DrugSess = 1;
Pat_Sub{7}.Fids = [ ] ; %NAS, LPA, RPA
Pat_Sub{7}.Dia = 'bvFTD';
Pat_Sub{7}.CarTask = 0;
Pat_Sub{7}.RSC = [0,0,0,0];
Pat_Sub{7}.EEG = [];

Pat_Sub{8}.Name = {'meg17_0012', 'meg17_0050', 'P8', 'DATA_0005.nii'};
Pat_Sub{8}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/21august/25585/20161121_U-ID37328/Series_005_t1_mprage_sag/';
Pat_Sub{8}.DrugSess = 1;
Pat_Sub{8}.Fids = [0.3 83.1 -27.0; -74.7 8.1 -56.6; 67.6 -7.6 -58.1 ] ; %NAS, LPA, RPA
Pat_Sub{8}.Dia = 'PSP';
Pat_Sub{8}.CarTask = 1;
Pat_Sub{8}.RSC = [1,0,1,1];
Pat_Sub{8}.EEG = [1,1];

Pat_Sub{9}.Name = {'meg17_0018', 'meg17_0034', 'P9', 'DATA_0002.nii'};
Pat_Sub{9}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/21august/25321/20151211_U-ID36619/Series_002_MPRAGE_1.25_iso/';
Pat_Sub{9}.DrugSess = 1;
Pat_Sub{9}.Fids = [3.6 92.0 -20.2; -76.3 5.9 -50.2; 72.8 10.7 -57.3] ; %NAS, LPA, RPA
Pat_Sub{9}.Dia = 'bvFTD';
Pat_Sub{9}.CarTask = 1;
Pat_Sub{9}.RSC = [1,1,1,1];
Pat_Sub{9}.EEG = [1,1];

Pat_Sub{10}.Name = {'meg17_0035', 'NA', 'P10', ''};
Pat_Sub{10}.Taldir = '24589 ';
Pat_Sub{10}.DrugSess = 2;
Pat_Sub{10}.Fids = [  ] ; %NAS, LPA, RPA
Pat_Sub{10}.Dia = 'bvFTD';
Pat_Sub{10}.CarTask = 0;
Pat_Sub{10}.RSC = [0,0,0,0];
Pat_Sub{10}.EEG = [];

Pat_Sub{11}.Name = {'meg17_0054', 'meg17_0074', 'P11', ''};
Pat_Sub{11}.Taldir = ''; %%no MRI?
Pat_Sub{11}.DrugSess = 2;
Pat_Sub{11}.Fids = [] ; %NAS, LPA, RPA
Pat_Sub{11}.Dia = 'bvFTD';
Pat_Sub{11}.CarTask = 0;
Pat_Sub{11}.RSC = [0,0,0,0];
Pat_Sub{10}.EEG = [];

Pat_Sub{12}.Name = {'NA', 'NA', 'P12', ''};
Pat_Sub{12}.Taldir = '';
Pat_Sub{12}.DrugSess = 0;
Pat_Sub{12}.Fids = [] ; %NAS, LPA, RPA
Pat_Sub{12}.Dia = '';
Pat_Sub{12}.CarTask = 0;
Pat_Sub{12}.RSC = [0,0,0,0]; % AD
Pat_Sub{10}.EEG = [];

Pat_Sub{13}.Name = {'meg17_0140', 'meg17_0152', 'P13', 'single_subj_T1.nii'};
Pat_Sub{13}.Taldir = 'tbc?'; %%
Pat_Sub{13}.DrugSess = 1;
Pat_Sub{13}.Fids = [2.2 83.8 -41.6;  -82.9 -18.9 -53.2; 85.6 -17.1 -53.2]; %NAS, LPA, RPA
Pat_Sub{13}.Dia = 'bvFTD';
Pat_Sub{13}.CarTask = 1;
Pat_Sub{13}.RSC = [1,0,1,0]; % s2 rsc2 too short
Pat_Sub{13}.EEG = [0,0];

Pat_Sub{14}.Name = {'meg17_0133', 'meg17_0135', 'P14', '25989_mprage_sag.nii'};
Pat_Sub{14}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/25989/20170714/mprage_sag/';
Pat_Sub{14}.DrugSess = 1;
Pat_Sub{14}.Fids = [7.0 83.9 -13.3; -72.0 6.5 -48.5; 80.2 -5.2 -48.5  ] ; %NAS, LPA, RPA
Pat_Sub{14}.Dia = 'bvFTD';
Pat_Sub{14}.CarTask = 1;
Pat_Sub{14}.RSC = [1,1,0,0];
Pat_Sub{14}.EEG = [1,0];

Pat_Sub{15}.Name = {'meg17_0158', 'meg17_0165', 'P15', 'DATA_0003.nii'};
Pat_Sub{15}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/21august/26565/20180207_U-ID39807/Series_003_Sag_MPRAGE/';
Pat_Sub{15}.DrugSess = 2;
Pat_Sub{15}.Fids = [-2.8 95.5 -14.9; -74.6 -15.6 -48.9; 75.3 -12.5 -52.0] ; %HOLLY NAS, LPA, RPA
Pat_Sub{15}.Dia = 'bvFTD';
Pat_Sub{15}.CarTask = 0;
Pat_Sub{15}.RSC = [0,0,1,1];
Pat_Sub{15}.EEG = [0,1];

Pat_Sub{16}.Name = {'NA', 'NA', 'P16', ''};
Pat_Sub{16}.Taldir = ''; %
Pat_Sub{16}.DrugSess = 0;
Pat_Sub{16}.Fids = [ ] ; %NAS, LPA, RPA
Pat_Sub{16}.Dia = 'PSP';
Pat_Sub{16}.CarTask = 0;
Pat_Sub{16}.RSC = [0,0,0,0];% Pulled out
Pat_Sub{16}.EEG = [];

Pat_Sub{17}.Name = {'meg17_177', 'meg17_0184', 'P17', '22118_mprage_125.nii'};
Pat_Sub{17}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/22118/';                                  
Pat_Sub{17}.DrugSess = 2;
Pat_Sub{17}.Fids = [ -1.6 99.3 -15.0; -74.3 -11.6 -42.2; 73.2 -4.0 -48.7] ; 
Pat_Sub{17}.Dia = 'bvFTD';
Pat_Sub{17}.CarTask = 0;
Pat_Sub{17}.RSC = [1,0,1,1];
Pat_Sub{17}.EEG = [0,0];

Pat_Sub{18}.Name = {'meg17_0192', 'meg17_0204', 'P18', 'DATA_0006.nii'};
Pat_Sub{18}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/25910/20171201_U-ID39328/Series_006_Sag_MPRAGE/';
Pat_Sub{18}.DrugSess = 2;
Pat_Sub{18}.Fids = [0.0 77.6 -8.1; -74.8 -9.4 -41.7; 76.0 -14.8 -47.1 ] ; %NAS, LPA, RPA
Pat_Sub{18}.Dia = 'bvFTD';
Pat_Sub{18}.CarTask = 1;
Pat_Sub{18}.RSC = [1,1,1,1];
Pat_Sub{18}.EEG = [1,1];

Pat_Sub{19}.Name = {'meg17_0219', 'meg17_0229', 'P19', 'DATA_0006.nii'};
Pat_Sub{19}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/p00363/26271/20171102_U-ID39110/Series_006_Sag_MPRAGE/';
Pat_Sub{19}.DrugSess = 1;
Pat_Sub{19}.Fids = [-0.6 80.0 -21.2; -80.0 -15.5 -52.1; 75.7 -7.4 -52.1  ] ; %NAS, LPA, RPA
Pat_Sub{19}.Dia = 'bvFTD';
Pat_Sub{19}.CarTask = 1;
Pat_Sub{19}.RSC = [1,1,1,1];
Pat_Sub{19}.EEG = [1,1];

Pat_Sub{20}.Name = {'NA', 'NA', 'P20', '.nii'};
Pat_Sub{20}.Taldir = '';
Pat_Sub{20}.DrugSess = 0;
Pat_Sub{20}.Fids = [] ; %NAS, LPA, RPA
Pat_Sub{20}.CarTask = 0;
Pat_Sub{20}.RSC = [0,0,0,0];
Pat_Sub{20}.EEG = [];

Pat_Sub{21}.Name = {'meg18_0006', 'meg18_0019', 'P21', 'DATA_0006.nii'};
Pat_Sub{21}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/21august/26461/20180117_U-ID39606/Series_006_Sag_MPRAGE/';
Pat_Sub{21}.DrugSess = 2;
Pat_Sub{21}.Fids = [0.0 79.3 -0.7; -72.1 -14.0 -40.9; 79.0 -1.2 -58.6 ] ; %NAS, LPA, RPA
Pat_Sub{21}.Dia = 'bvFTD';
Pat_Sub{21}.CarTask = 1;
Pat_Sub{21}.RSC = [1,1,1,1];
Pat_Sub{21}.EEG = [1,1];

Pat_Sub{22}.Name = {'meg18_0055', 'meg18_0058', 'P22', 'DATA_0008.nii'};
Pat_Sub{22}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/21august/27216/20180704_U-ID41303/Series_008_mp2rage_sag_p3_0.75mm_UNI_Images/';
Pat_Sub{22}.DrugSess = 1;
Pat_Sub{22}.Fids = [2.1 77.8 -13.8; -75.7 -9.6 -41.7; 71.4 -3.6 -50.2] ; %NAS, LPA, RPA SCAN ODD
Pat_Sub{22}.Dia = 'bvFTD';
Pat_Sub{22}.CarTask = 1;
Pat_Sub{22}.RSC = [1,1,1,1];
Pat_Sub{22}.EEG = [1,1];

Pat_Sub{23}.Name = {'NA', 'NA', 'P23', ''};
Pat_Sub{23}.Taldir = '';
Pat_Sub{23}.DrugSess = 0;
Pat_Sub{23}.Fids = [] ; %NAS, LPA, RPA
Pat_Sub{23}.Dia = 'PSP';
Pat_Sub{23}.CarTask = 0;
Pat_Sub{23}.RSC = [0,0,0,0];
Pat_Sub{23}.EEG = [];

Pat_Sub{24}.Name = {'meg18_0067', 'meg18_0071', 'P24', 'DATA_0005.nii'};
Pat_Sub{24}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/10oct/26418/20171218_U-ID39485/Series_005_Sag_MPRAGE_ND/';
Pat_Sub{24}.DrugSess = 1;
Pat_Sub{24}.Fids = [-1.3 84.6 -8.9; -76.3 -5.1 -49.5; 78.7 5.7 -57.2 ] ; %NAS, LPA, RPA
Pat_Sub{24}.Dia = 'PSP';
Pat_Sub{24}.CarTask = 1;
Pat_Sub{24}.RSC = [1,0,1,1]; % s1 rsc2 not doing MF trans for unknown reason, come back to at later date
Pat_Sub{24}.EEG = [1,0];

Pat_Sub{25}.Name = {'meg18_0068', 'meg18_0076', 'P25', 'DATA_0006.nii'};
Pat_Sub{25}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/21august/26907/20180424_U-ID40557/Series_006_Sag_MPRAGE/';
Pat_Sub{25}.DrugSess = 2;
Pat_Sub{25}.Fids = [-8.0 86.9 0.0; -86.5 13.8 -51.8; 71.8 12.4 -55.8] ; %NAS, LPA, RPA
Pat_Sub{25}.Dia = 'bvFTD';
Pat_Sub{25}.CarTask = 1;
Pat_Sub{25}.RSC = [1,1,1,1];
Pat_Sub{25}.EEG = [1,1];

Pat_Sub{26}.Name = {'meg18_0078', 'meg18_0084', 'P26', 'DATA_0008.nii'};
Pat_Sub{26}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/21august/25799/20180817_U-ID41734/Series_008_mp2rage_sag_p3_0.75mm_UNI_Images/';
Pat_Sub{26}.DrugSess = 2;
Pat_Sub{26}.Fids = [3.6 79.8 -5.9; -84.6 -18.4 -37.7; 72.2 -13.5 -45.1 ] ; %ODD SCAN %NAS, LPA, RPA
Pat_Sub{26}.Dia = 'PSP';
Pat_Sub{26}.CarTask = 1;
Pat_Sub{26}.RSC = [1,1,1,1];
Pat_Sub{26}.EEG = [1,1];

Pat_Sub{27}.Name = {'NA', 'NA', 'P27', 'DATA_0008.nii'};
Pat_Sub{27}.Taldir = '\imaging\na01\HPHI_MRIs\tallie\10oct\27389\20180816_U-ID41705\Series_008_mp2rage_sag_p3_0.75mm_UNI_Images\';
Pat_Sub{27}.DrugSess = 2;
Pat_Sub{27}.Fids = [] ; %NAS, LPA, RPA
Pat_Sub{27}.Dia = 'PSP';
Pat_Sub{27}.CarTask = 0;
Pat_Sub{27}.RSC = [0,0,0,0];
Pat_Sub{27}.EEG = [];

Pat_Sub{28}.Name = {'meg18_0116', 'meg18_0123', 'P28', 'single_subj_T1.nii'};
Pat_Sub{28}.Taldir = '';
Pat_Sub{28}.DrugSess = 1;
Pat_Sub{28}.Fids = [2.2 83.8 -41.6;  -82.9 -18.9 -53.2; 85.6 -17.1 -53.2] ; %NAS, LPA, RPA
Pat_Sub{28}.Dia = 'bvFTD';
Pat_Sub{28}.CarTask = 1;
Pat_Sub{28}.RSC = [1,0,0,0];
Pat_Sub{28}.EEG = [1,0];

Pat_Sub{29}.Name = {'meg18_0118', 'meg18_0129', 'P29', 'DATA_0008.nii'};
Pat_Sub{29}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/10oct/26413/20180801_U-ID41545/Series_008_mp2rage_sag_p3_0.75mm_INV2/';
Pat_Sub{29}.DrugSess = 2;
Pat_Sub{29}.Fids = [3.4 77.5 -9.9; -72.7 -7.3 -44.8; 67.9 -9.4 -48.0] ; %NAS, LPA, RPA
Pat_Sub{29}.Dia = 'PSP';
Pat_Sub{29}.CarTask = 1;
Pat_Sub{29}.RSC = [1,0,1,0];
Pat_Sub{29}.EEG = [1,1];

Pat_Sub{30}.Name = {'meg19_0024', 'meg19_0050', 'P30', '27449_mp2rage.nii'};
Pat_Sub{30}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/';
Pat_Sub{30}.DrugSess = 2;
Pat_Sub{30}.Fids = [3.6 135.2 -13.1; -68.0 40.0 -58.7; 73.9 51.7 -62.6] ; %NAS, LPA, RPA
Pat_Sub{30}.Dia = 'PSP';
Pat_Sub{30}.CarTask = 1;
Pat_Sub{30}.RSC= [1,1,1,1];
Pat_Sub{30}.EEG = [1,1];

Pat_Sub{31}.Name = {'meg18_0127', 'meg18_0139', 'P31', 'DATA_0008.nii'};
Pat_Sub{31}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/10oct/26135/20180913_U-ID41900/Series_008_mp2rage_sag_p3_0.75mm_UNI_Images/';
Pat_Sub{31}.DrugSess = 2;
Pat_Sub{31}.Fids = [2.6 81.8 13.8; -75.6 -1.3 -39.3; 65.1 -13.0 -35.4 ] ; %NAS, LPA, RPA
Pat_Sub{31}.Dia = 'PSP';
Pat_Sub{31}.CarTask = 1; %% 
Pat_Sub{31}.RSC= [1,0,1,1];
Pat_Sub{31}.EEG = [1,1];

Pat_Sub{32}.Name = {'meg18_0135', 'meg18_0146', 'P32', '25935_mp2rage.nii'};
Pat_Sub{32}.Taldir = '/imaging/na01/HPHI_MRIs/11Jan19/';
Pat_Sub{32}.DrugSess = 1;
Pat_Sub{32}.Fids = [-2.3 82.1 -8.7; -81.7 -24.5 -43.7; 73.8 -20.0 -52.8] ; %NAS, LPA, RPA
Pat_Sub{32}.Dia = 'PSP';
Pat_Sub{32}.CarTask = 1; %% will be a 1 if completed
Pat_Sub{32}.RSC= [1,1,1,1];
Pat_Sub{32}.EEG = [1,0];

Pat_Sub{33}.Name = {'meg18_0143', 'NA', 'P33', '25550_mp2rage.nii'}; %sick on both days
Pat_Sub{33}.Taldir = '/imaging/na01/HPHI_MRIs/11Jan19/';
Pat_Sub{33}.DrugSess = 1;
Pat_Sub{33}.Fids = [6.7 131.2 -11.7; -71.2 27.7 -50.6; 78.6 19.5 -43.6] ; % HOLLY NAS, LPA, RPA
Pat_Sub{33}.Dia = 'PSP';
Pat_Sub{33}.CarTask = 0; %% %sick on both days
Pat_Sub{33}.RSC= [1,0,0,0];
Pat_Sub{33}.EEG = [1,0];

Pat_Sub{34}.Name = {'meg18_0147', 'meg18_0163', 'P34', '27626_mp2rage.nii'};
Pat_Sub{34}.Taldir = '/imaging/na01/HPHI_MRIs/11Jan19/';
Pat_Sub{34}.DrugSess = 1;
Pat_Sub{34}.Fids = [-3.8 68.5 -28.0; -75.9 -30.7 -30.6; 64.4 -29.4 -33.2] ; %NAS, LPA, RPA
Pat_Sub{34}.Dia = 'PSP';
Pat_Sub{34}.CarTask = 1; %% will be a 1 if completed
Pat_Sub{34}.RSC= [1,1,1,1];
Pat_Sub{34}.EEG = [1,1];

Pat_Sub{35}.Name = {'meg18_0155', 'meg18_0167', 'P35', '26431_mp2rage.nii'};
Pat_Sub{35}.Taldir = '/imaging/na01/HPHI_MRIs/11Jan19/';
Pat_Sub{35}.DrugSess = 2;
Pat_Sub{35}.Fids = [22.8 118.3 -6.3;-71.1 24.3 -12.6; 71.6 18.0 -13.9] ; % HOLLY NAS, LPA, RPA
Pat_Sub{35}.Dia = 'bvFTD';
Pat_Sub{35}.CarTask = 0; %% not able to complete AT
Pat_Sub{35}.RSC= [1,1,1,0];
Pat_Sub{35}.EEG = [1,0];

Pat_Sub{36}.Name = {'meg18_0165', 'NA', 'P36', '26949_mp2rage.nii'};
Pat_Sub{36}.Taldir = '/imaging/na01/HPHI_MRIs/11Jan19/';
Pat_Sub{36}.DrugSess = 1;
Pat_Sub{36}.Fids = [-12.6 119.7 -18.2; -65.0 5.4 -68.1; 70.8 13.7 -40.8] ; %HOLLY NAS, LPA, RPA
Pat_Sub{36}.Dia = 'PSP';
Pat_Sub{36}.CarTask = 0; %% sick on day not able to complete
Pat_Sub{36}.RSC= [1,1,0,0];
Pat_Sub{36}.EEG = [1,0];

Pat_Sub{37}.Name = {'meg19_0015', 'meg19_0038', 'P37', '26796_mp2rage.nii'};
Pat_Sub{37}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/Feb19/';
Pat_Sub{37}.DrugSess = 2;
Pat_Sub{37}.Fids = [1.1 68.8 -32.0; -70.3 -29.1 -35.5 ; 67.9 -21.1 -44.7] ; %NAS, LPA, RPA
Pat_Sub{37}.Dia = 'PSP';
Pat_Sub{37}.CarTask = 1; 
Pat_Sub{37}.RSC= [1,1,1,1];
Pat_Sub{37}.EEG = [1,1];

Pat_Sub{38}.Name = {'meg19_0042', 'meg19_0069', 'P38', '25382_20181206.nii'};
Pat_Sub{38}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/Feb19/';
Pat_Sub{38}.DrugSess = 1;
Pat_Sub{38}.Fids = [15.1 130.6 4.6; -76.1 35 -28; 69.5 18.6 -53.9] ; %NAS, LPA, RPA
Pat_Sub{38}.Dia = 'PSP';
Pat_Sub{38}.CarTask = 0; %% no second session for AT
Pat_Sub{38}.RSC= [1,1,1,0];
Pat_Sub{38}.EEG = [1,1];

Pat_Sub{39}.Name = {'meg19_0062', 'meg19_0077', 'P39', '27928_20190129_mprage.nii'};
Pat_Sub{39}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/Feb19/';
Pat_Sub{39}.DrugSess = 2;
Pat_Sub{39}.Fids = [4.9 78.8 -5.2; -72.6 -22.1 -61.1 ; 73.5 -25.9 -58.6] ; %NAS, LPA, RPA 3T
Pat_Sub{39}.Dia = 'bvFTD';
Pat_Sub{39}.CarTask = 1; %% very good at AT
Pat_Sub{39}.RSC= [1,1,1,1];
Pat_Sub{39}.EEG = [1,1];

Pat_Sub{40}.Name = {'meg19_0115', 'meg19_0136', 'P40', '27928_20190129_mprage.nii'};
Pat_Sub{40}.Taldir = '/imaging/na01/HPHI_MRIs/tallie/Feb19/';
Pat_Sub{40}.DrugSess = 2; %%%% CHANGE WHEN I KNOW!!!
Pat_Sub{40}.Fids = [-0.6 79.3 -10; -78.8 -26.1 -55.7; 72.4 -40.1 -55.7] ; % HOLLY NAS, LPA, RPA
Pat_Sub{40}.Dia = 'bvFTD';
Pat_Sub{40}.CarTask = 0; %% 
Pat_Sub{40}.RSC= [1,1,1,0]; % Error on MF trans for s2 rsc2, not sure why atm
Pat_Sub{40}.EEG = [1,1];


%% addin MRI dir and MEG dir for AT
for n = 1:length(Pat_Sub)
    Pat_Sub{n}.MRIdir = sprintf('/imaging/lh01/FTD2015/MRIs/%s/',Pat_Sub{n}.Name{3} );
    Pat_Sub{n}.megpath{1} = sprintf('%s%s/MF/SPM12analysis/', wd, Pat_Sub{n}.Name{1});
    Pat_Sub{n}.megpath{2} = sprintf('%s%s/MF/SPM12analysis/', wd, Pat_Sub{n}.Name{2});
end


%% use this to make dir of MRIs
maknewdirs = 0;
if maknewdirs
    includingsubs = [37 39];
    %waitingfor = [13 24];
    
    for nnn = 1:length(includingsubs) %subnum
        n= includingsubs(nnn);%which subs to run
        filename = char(Pat_Sub{n}.Name(4));
        subname = char(Pat_Sub{n}.Name(3));
        taldir = Pat_Sub{n}.Taldir;
        
        mkdir_text = sprintf('mkdir ''%s''',subname);
        cd_text_tal = sprintf('cd ''%s''', taldir);
        copyfile_text = sprintf('copyfile (''%s'', ''/imaging/lh01/FTD2015/MRIs/%s/'')', filename, subname);
        
        cd '/imaging/lh01/FTD2015/MRIs'
        eval(mkdir_text)
        eval(cd_text_tal)
        eval(copyfile_text)
        
    end
    
end







end %% end function