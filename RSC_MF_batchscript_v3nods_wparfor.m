% RSC Max Filter 2016 LEH - adapted by HNP 2019
% v2 - standard trans to default space, v3 - trans to 1st sess then to
% default
% use this 
clear all
addpath /neuro/meg_pd_1.2
addpath /imaging/local/meg_misc/
addpath(genpath('/imaging/hp02/spm12b'));
addpath /imaging/hp02/TGB/matlab_scripts

cd /imaging/hp02/TGB/rest_closed/

%% Matlab pool stuff
 % close any cbupool currently open
        cbupool(40); %11 or 5

%% Set up groups

%% Controls
con_data = Control_Subs_Data;
for ss = 1:length(con_data)
    DrugSess = con_data{1,ss}.DrugSess;
    PlcbSess = 3-DrugSess;
    cons_drug{ss} = con_data{1,ss}.Name{DrugSess};
    cons_plcb{ss} = con_data{1,ss}.Name{PlcbSess};
end

cons_MF = [cons_drug,cons_plcb];

% Cons_s1={ 'meg16_0250', 'meg17_0036', 'meg17_0068',	'meg17_0067', 'meg17_0020', ...
%     'meg17_0048', 'meg17_0084', 'meg17_0089', 'meg17_0098', 'meg17_0101', ...
%     'meg17_0124', 'meg17_0118', 'meg17_0113', 'meg17_0123', 'meg17_0126', ...
%     'meg17_0179',	'meg17_0198', 'meg17_0208', 'meg17_215', 'meg17_0052' };
% Cons_s2={ 'meg16_0261', 'meg17_0016', 'meg17_0038', 'meg17_0030', 'meg17_0044', ...
%     'meg17_0091', 'meg17_0057', 'meg17_0115', 'meg17_0081', 'meg17_0096', ...
%     'meg17_0119', 'meg17_0122', 'meg17_0156', 'meg17_0114', 'meg17_0147', ...
%     'meg17_0189', 'meg17_0210', 'meg17_0195', 'meg17_0231', 'meg17_0077'};
%% Patients
pat_data = Patient_Subs_Data;
for ss = 1:length(pat_data)
    DrugSess = pat_data{1,ss}.DrugSess;
    PlcbSess = 3-DrugSess;
    if DrugSess ~= 0
    pats_drug{ss} = pat_data{1,ss}.Name{DrugSess};
    pats_plcb{ss} = pat_data{1,ss}.Name{PlcbSess};
    else
        pats_drug{ss} = 'NA';
        pats_plcb{ss} = 'NA';
    end
end

% Select those to be analysed
j = 1;ss = 1;
while ss <= 40
    if ~strcmp(pats_drug{ss}, 'NA')
        pats_MF{j} = pats_drug{ss};
        j = j+1;
    end
    
    if ~strcmp(pats_plcb{ss}, 'NA')
        pats_MF{j} = pats_plcb{ss};
        j = j+1;
    end
    
    ss = ss+1;
end

% Pats = {'meg16_0284',  'meg16_0295',  'meg16_0287',  'meg16_0308',  'meg16_0335', ... 	  
%         'NA',          'meg17_0002',  'meg17_0060',  'meg17_0005',  'meg17_0028',  ...
%        'meg17_0007',   'meg17_0039',  'meg17_0012',  'meg17_0050',  'meg17_0018',  ...
%        'meg17_0034',   'meg17_0035',  'NA',          'meg17_0054',  'meg17_0074',  ...
%        'NA',           'NA',          'meg17_0140',  'meg17_0152',  'meg17_0133',  ...
%        'meg17_0135',   'meg17_0158',  'meg17_0165',  'NA',          'NA',      ... 
%        'meg17_177',    'meg17_0184',  'meg17_0192',  'meg17_0204',  'meg17_0219',  ...
%        'meg17_0229',   'NA',          'NA',          'meg18_0006',  'meg18_0019',  ...
%        'meg18_0055',   'meg18_0058',  'NA',          'NA',          'meg18_0067',  ...
%        'meg18_0071',   'meg18_0068',  'meg18_0076',  'meg18_0078',  'meg18_0084', ...
%        'meg18_0116',   'meg18_0129',  'meg18_0123',  'meg18_0118',  'meg18_0127', ...
%        'meg18_0139', 'meg18_0135', 'meg18_0146', 'meg18_0147', 'meg18_0155'}; 
%    
% NewPats = {'meg18_0076', 'meg19_0024', 'meg19_0050', 'meg18_0127', 'meg18_0139', 'meg18_0135', 'meg18_0146', ...
% 'meg18_0147', 'meg18_0163', 'meg19_0015', 'meg19_0038', 'meg19_0062', 'meg19_0077'};
%p2:p26



%pats_to_include = 1 %1:12 %[1:4, 9:10, 13:16, 23:26, 33:36, 39:42, 45:50  51:55] ;
%singles = [7, 17 , 32 ]; %consider doing later

% meg_badchans_cs1 = { ['0813 ', '2113 ', '2122 ']; ['0813 ', ];  ['0813 ']}    ;
% meg_badchans_cs2 = { ['0813 ']; ['0813 ', ];  ['0813 ']}    ;
% meg_badchans_ps1 = { ['0813 ']; ['0813 ', ];  ['0813 ']}    ;
% meg_badchans_ps2 = { ['0813 ']; ['0813 ', ];  ['0813 ']}    ;

subs2proc = cons_MF;  %Pats % [Cons_s1 Cons_s2];
% meg_badchans = ['0813 '] %[meg_badchans_cs1 meg_badchans_cs2];

dosubs =1:length(cons_MF) ;

%% Set UP
%which MF process to run
badchanflag = 1;
maxsssflag = 1;
maxtransflag = 0;
chancheckflag = 0;

types{1} = '1';
types{2} = '2';
%types{3} = '3';
%types{4} = '4'; %if 4th need to re look at data sets

%settings
trans_offset_maxf = [0 -13 +6];

clear cmdb; clear cmd_badgo; clear cmd; clear badchstr; clear badch
clear filename; clear headpoints; clear fit
%% Maxfilter
disp('Bad Channels')
for ss = 1:length(dosubs)
    
    n = dosubs(ss)
    
    %make a dir to put the new MF sss files in
    megpath = strcat('/imaging/hp02/TGB/rest_closed/MF/Controls/', subs2proc{n},'/');
    if exist(megpath)~=7; mkdir(megpath);end
    
    % *** raw data in cbu dir **
    rp1 = strcat('/megdata/cbu/ftdrug/',subs2proc{n}, '/'); %
    rp2 = dir(rp1);
    rawpath = strcat(rp1, rp2(3).name, '/');
 
    allraw = spm_select('List', rawpath, '.*raw.fif'); %% all files
    allraw1 = lower(allraw); %change to lowe case
    allraw_rsc = [];
    for nnn = 1:size(allraw1,1)
        if strfind(allraw1(nnn,:),'rsc')
            allraw_rsc = [allraw_rsc; allraw(nnn,:)];
        end
    end
    allraw = allraw_rsc;
    Nses{ss} = 2;%size(allraw,1); %number of files
    

     %make a dir to put the new MF sss files in
     cd (megpath)
    rfxrootdir = fullfile(megpath);
%     if exist(rfxrootdir)~=7; mkdir(megpath,'/MF/');end
%     cd('MF')
    
    %************************************************************%
    %                 Loop for sessions                          %
    %************************************************************%
    for k = 1:Nses{ss} %raw file to process (sess1 etc)
        
        
        clear cmdb;  clear cmd; clear badchstr; clear badch
        clear headpoints; clear fit
        
        % overwrite file names so that they include the mistyped names
        load('/imaging/hp02/TGB/rest_closed/MF/control_filenames_inclerrornames.mat');
        
        %filename{ss}.(sprintf('k%s',num2str(k))) = strcat(rawpath, allraw(k,:) ); %1st raw file in dir
        currfile{ss}.(sprintf('k%s',num2str(k))) = types{k};
        
        
        %************************************************************%
        % Danny Mitchell & JT: headpos code Fit sphere to HPI points %
        %************************************************************%
        
        %*** option for head points from riks pipeline ***
        %
        headptsfile  = sprintf('%s%s_headpoints_%s.txt',rfxrootdir,subs2proc{n},currfile{ss}.(sprintf('k%s',num2str(k))));
        [co ki] = hpipoints(filename{ss}.(sprintf('k%s',num2str(k))));
        headpoints = co(:,ki>1)'; % don't include the fiducial points
        headpoints = headpoints(~(headpoints(:,2)>0 & headpoints(:,3)<0),:);% Remove nose points
        cd(megpath)
        save('hpipoints.txt', '-ASCII', 'headpoints')
        cmd_fit = '/neuro/bin/util/fit_sphere_to_points hpipoints.txt';% Fit sphere:
        [status spherefit] = unix(cmd_fit);
        if length(spherefit) < 1;
            error('Spherefit failed!')
        end
        fit = str2num(spherefit)*1000; % m to mm;
        fit1  = num2str(fit(1));
        fit2  = num2str(fit(2));
        fit3  = num2str(fit(3));
        % use in -frame head -origin
        allfits{ss}.(sprintf('k%s',num2str(k)))  = [fit1, ' ', fit2, ' ', fit3] ;
        
        %************************************************************%
        %                 get  bad channels                          %
        %************************************************************%
        
        if  badchanflag
            cmdb{1} =sprintf(' /neuro/bin/util/maxfilter-2.2 ');
            cmdb{2} =sprintf(' -f %s',filename{ss}.(sprintf('k%s',num2str(k))) );
            cmdb{3} =sprintf(' -o %s/%s_rsc%s_badCH.fif ',megpath,subs2proc{n},currfile{ss}.(sprintf('k%s',num2str(k))));
            cmdb{4} =sprintf(' -ctc /neuro/databases/ctc/ct_sparse.fif ');
            cmdb{5} =sprintf(' -cal /neuro/databases/sss/sss_cal.dat ');
            %cmdb{6} =sprintf(' -ds 2 -format short -hpistep 1000 -hpisubt amp -in 8  -out 3 ');
            cmdb{7} =sprintf(' -force -frame head -origin ');
            cmdb{8} = sprintf(' %s %s %s', fit1, fit2, fit3);
            cmdb{9} = sprintf(' -autobad 20' ); %-skip 20.0 9999.0 ');
            cmdb{10} =sprintf(' -v | tee %s%s_rsc%s_bad_logfile.log ', megpath,subs2proc{n},currfile{ss}.(sprintf('k%s',num2str(k))));
            
            cmd_badgo{ss}.(sprintf('k%s',num2str(k))) = strcat(cmdb{1}, cmdb{2}, cmdb{3}, cmdb{4}, cmdb{5}, cmdb{7}, cmdb{8}, cmdb{9}, cmdb{10} );
            
        end
    end %ses
end %subs

%% Bad channels
if badchanflag
parfor ss = 1:length(dosubs)
    for k = 1:Nses{ss}
        
        unix(cmd_badgo{ss}.(sprintf('k%s',num2str(k))));
    end
end
end

%% SSS and Movecomp

disp('SSS and Movecomp')
clear sss1donealready sss2donealready

for ss = 1:length(dosubs)
    for k = 1:Nses{ss} %raw file to process (sess1 etc)
        n = dosubs(ss)
        
        %
        megpath = strcat('/imaging/hp02/TGB/rest_closed/MF/Controls/', subs2proc{n},'/');
        megpath_pf{ss}.(sprintf('k%s',num2str(k))) = strcat('/imaging/hp02/TGB/rest_closed/MF/Controls/', subs2proc{n},'/');
        
        %get bad channels into a variable as strings from log file and save as
        %sub_badchanfile
        getbadch= sprintf('cat %s%s_rsc%s_bad_logfile.log | sed -n ''''/Static/p'''' | cut -f 5- -d '' '' | tee %s%s_rsc%s_badchanfile.txt ', megpath,subs2proc{n},currfile{ss}.(sprintf('k%s',num2str(k))),megpath,subs2proc{n},currfile{ss}.(sprintf('k%s',num2str(k))));
        [dnr, bcs] = unix(getbadch);
        allbadch=str2num(bcs); %turn into nums
        
        badchstr=['0813 ' ]; % [] to include listed badchans from recording
        for xyz =1:size(allbadch,2) %individual variables into strings
            badch{xyz} = num2str(allbadch(1,xyz));
            badchstr = [badchstr, ' ', badch{xyz}, ' '];
            badchstr_save{ss}.bad = badchstr;
            badchstr_save{ss,k}.file = subs2proc{n}; 
        end

        %% meg16_0287
       % badchstr=[ '631 1711 1721 1831 1841 1931 2021 2041 2131 2311 2331 2541' ];
        %************************************************************%
        %                      MF SSS                                %
        %************************************************************%
        %trying to automate sss with movecomp and st.....
        %may have hpi failure so do without movecomp or add inter to interpolate
        %head pos if failures are few and short.
        
        % load corrected filenames:
        load('/imaging/hp02/TGB/rest_closed/MF/control_filenames_inclerrornames.mat');
        
        inputfile =  filename{ss}.(sprintf('k%s',num2str(k))) ; %sss 1 with movecomp
        sssfile1{ss}.(sprintf('k%s',num2str(k))) = sprintf('rsc%s_sss_%s_wods.fif',currfile{ss}.(sprintf('k%s',num2str(k))), subs2proc{n});
        %adds ST after movecomp
        sssfile2{ss}.(sprintf('k%s',num2str(k))) = sprintf('rsc%s_sss_%s.fif',currfile{ss}.(sprintf('k%s',num2str(k))), subs2proc{n});
      
        
        if maxsssflag
            cmd{1} ='/neuro/bin/util/maxfilter-2.2 ';
            cmd{2} =sprintf(' -f %s', filename{ss}.(sprintf('k%s',num2str(k))));
            cmd{3} =sprintf(' -o %s/%s', megpath, sssfile1{ss}.(sprintf('k%s',num2str(k))));
            cmd{4} =sprintf(' -hp %s/rsc_sss_new_%s_%s.pos', megpath, subs2proc{n}, currfile{ss}.(sprintf('k%s',num2str(k))));
            cmd{5} = sprintf (' -autobad off -bad %s ', badchstr);%badchstr ; % for SSS
            cmd{6} =' -force -frame head -origin ';
            cmd{7} = sprintf(' %s', allfits{ss}.(sprintf('k%s',num2str(k))) );
            %cmd{7} = sprintf(' %s %s %s', fit1, fit2, fit3);
            cmd{8} = ' -ctc /neuro/databases/ctc/ct_sparse.fif ';
            cmd{9} = ' -cal /neuro/databases/sss/sss_cal.dat ';
            cmd{10} = ' -st  ' ;  %-skip 0 20.0 -hpisubt amp -hpistep 200-format float%hpisubt removed high freq hpi
            cmd{11} = '  -movecomp inter  '; %-movecomp %inter %-skip 0 20.0 - not working
            cmd{12} =sprintf(' -v | tee %s_%s_ssslogfile.log ', subs2proc{n}, currfile{ss}.(sprintf('k%s',num2str(k))) );
            
            %movecomp no ST
            cmd_sssgo{ss}.(sprintf('k%s',num2str(k))) = strcat(cmd{1}, cmd{2}, cmd{3}, cmd{5}, cmd{6}, cmd{7}, cmd{8}, cmd{9}, cmd{11}, cmd{12});
            % with ST (based on movecomp file)
            cmd{13} =sprintf(' -f %s/%s', megpath, sssfile1{ss}.(sprintf('k%s',num2str(k))));
            cmd{14} =sprintf(' -o %s/%s', megpath, sssfile2{ss}.(sprintf('k%s',num2str(k))));
            cmd_sssgo_wST{ss}.(sprintf('k%s',num2str(k))) = strcat(cmd{1}, cmd{13}, cmd{14}, cmd{5}, cmd{6}, cmd{7}, cmd{8}, cmd{9}, cmd{10}, cmd{12});
           
            %if movecomp fails with ST
            cmd_sssgo_nmc{ss}.(sprintf('k%s',num2str(k))) = strcat(cmd{1}, cmd{2}, cmd{14}, cmd{5}, cmd{6}, cmd{7}, cmd{8}, cmd{9},cmd{10}, cmd{12});
         
            %unix(cmd_sssgo)
            cd(sprintf('%s/',megpath))
            sss1donealready{ss}.(sprintf('k%s',num2str(k))) = numel(dir(sssfile1{ss}.(sprintf('k%s',num2str(k))) ) );
            sss2donealready{ss}.(sprintf('k%s',num2str(k))) = numel(dir(sssfile2{ss}.(sprintf('k%s',num2str(k))) ) );
        end
                
    end %ses
end %subs

cd /imaging/hp02/TGB/rest_closed/MF/Controls
%save('badChans_controls.mat', 'badchstr_save');
%%
if maxsssflag
    %% SSS w movecomp
    parfor ss = 1:length(dosubs)
        ss
       
        for k = 1:Nses{ss}
            if sss1donealready{ss}.(sprintf('k%s',num2str(k))) < 1
                cd(megpath_pf{ss}.(sprintf('k%s',num2str(k))) )
                unix(cmd_sssgo{ss}.(sprintf('k%s',num2str(k))));
                % create movement figure
                %movement_figure(subs2proc{ss},k);
                %close all
            end
        end
    end
 %  
    % sss with ST
    parfor ss = 1:length(dosubs)
        for k = 1:Nses{ss}
            if sss2donealready{ss}.(sprintf('k%s',num2str(k))) < 1
                cd(megpath_pf{ss}.(sprintf('k%s',num2str(k))) )
                %cd('MF')
                unix(cmd_sssgo_wST{ss}.(sprintf('k%s',num2str(k))));
                
                
            end
        end
    end
   
    % if movecomp failed try with ST
    fileswithoutmc = [];
    for ss = 1:length(dosubs)
        ss
        for k = 1:Nses{ss}
            cd(megpath_pf{ss}.(sprintf('k%s',num2str(k))) )
            %cd('MF')
            if exist( sssfile1{ss}.(sprintf('k%s',num2str(k))) )
                
            else
                unix(cmd_sssgo_nmc{ss}.(sprintf('k%s',num2str(k))));
                fileswithoutmc = [fileswithoutmc; ss, k];
            end
            
        end
    end
  end
% Subject c17 s1 rsc2 and s2 rsc2 have too many bad channels and will not
% maxfilter atm. So need to set Nses to 1:
Nses{17} = 1; Nses{37} = 1;
%% Now Trans to 2nd fif file then Default 
%take  output file (*sss*.fif) and trans it.

if maxtransflag
    for ss = 1:11%length(dosubs)
        for k = 1:Nses{ss} %raw file to process (sess1 etc)
            if Nses{ss}>1 && k>1%%trans to first session if it exists
                n = dosubs(ss)

                megpath = strcat('/imaging/hp02/TGB/rest_closed/MF/Controls/', subs2proc{n},'/');
                              
                inputfile =  sssfile2{ss}.(sprintf('k%s',num2str(k)));
                transfile{ss}.(sprintf('k%s',num2str(k))) = sprintf('rsc%s_transA_%s.fif',currfile{ss}.(sprintf('k%s',num2str(k))), subs2proc{n});
                transtofile = sssfile1{ss}.k1; %trans to first session if it exists
                
                cmd{1} =sprintf(' /neuro/bin/util/maxfilter-2.2 ');
                cmd{2} =sprintf(' -f  %s/%s ',megpath, inputfile );
                cmd{3} =sprintf(' -o %s/%s ',megpath, transfile{ss}.(sprintf('k%s',num2str(k))));
                cmd{4} =sprintf(' -trans  %s ', transtofile );
                %cmd{5} = sprintf ('-autobad off -bad 0612 %s ', badchstr);%badchstr ; % for SSS
                cmd{6} =' -force -frame head -origin ';
                cmd{7} = sprintf(' %s %s %s', fit1, fit2, fit3);
                %cmd{5} =sprintf(' -ctc /neuro/databases/ctc/ct_sparse.fif ');
                %cmd{6} =sprintf(' -cal /neuro/databases/sss/sss_cal.dat ');
                %cmd{7} = sprintf(' -skip 0 20.0 ');
                cmd{8} =sprintf(' -v | tee %s/rsc_logfile_trans_%s_%s.log ', megpath, subs2proc{n},currfile{ss}.(sprintf('k%s',num2str(k))));
                
                cmd_transgo{ss}.(sprintf('k%s',num2str(k))) = strcat(cmd{1}, cmd{2}, cmd{3}, cmd{4}, cmd{6}, cmd{7}, cmd{8});
                %unix(cmd_transgo)
                
            end %nses
        end %if > 1session
        
        
    end %subs
    %% trans to default
    for ss = 1:11%length(dosubs)
        for k = 1:Nses{ss} %raw file to process (sess1 etc)
            n = dosubs(ss)

            megpath = strcat('/imaging/hp02/TGB/rest_closed/MF/Controls/', subs2proc{n},'/');
            if k>1 && Nses{ss} > 1
                inputfile =  transfile{ss}.(sprintf('k%s',num2str(k)));
            else
                inputfile =  sssfile2{ss}.(sprintf('k%s',num2str(k)));
            end
            transfile2 = sprintf('rsc%s_trans_%s.fif',currfile{ss}.(sprintf('k%s',num2str(k))), subs2proc{n});
                        
            cmd{1} =sprintf(' /neuro/bin/util/maxfilter-2.2 ');
            cmd{2} =sprintf(' -f  %s/%s ',megpath, inputfile );
            cmd{3} =sprintf(' -o %s/%s ',megpath, transfile2);
            cmd{4} =sprintf(' -trans default ');
            %cmd{5} = sprintf ('-autobad off -bad 0612 %s ', badchstr);%badchstr ; % for SSS
            cmd{6} =' -force -frame head -origin ';
            cmd{7} = sprintf(' %s %s %s', fit1, fit2, fit3);
            %cmd{5} =sprintf(' -ctc /neuro/databases/ctc/ct_sparse.fif ');
            %cmd{6} =sprintf(' -cal /neuro/databases/sss/sss_cal.dat ');
            %cmd{7} = sprintf(' -skip 0 20.0 ');
            cmd{8} =sprintf(' -v | tee %s/rsc_logfile_trans_%s_%s.log ', megpath, subs2proc{n},currfile{ss}.(sprintf('k%s',num2str(k))));
            
            cmd_transgo2{ss}.(sprintf('k%s',num2str(k))) = strcat(cmd{1}, cmd{2}, cmd{3}, cmd{4}, cmd{6}, cmd{7}, cmd{8});
            %unix(cmd_transgo)
            
        end
        
        
        
    end %ses
    
    
end %subs
%%
for ss = 11%:length(dosubs)
    for k = 1:Nses{ss}
             cd(megpath_pf{ss}.(sprintf('k%s',num2str(k))) )
        %cd('MF')
        if k >1
              unix(cmd_transgo{ss}.(sprintf('k%s',num2str(k))) );
        end
        unix(cmd_transgo2{ss}.(sprintf('k%s',num2str(k))) )
    end
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% Olaf's channel checking (have not used for a while; should put in SVN!)
% if chancheckflag
%     addpath /imaging/olaf/MEG/artefact_scan_tool/
%     global fiffs tasks values criteria sensors datapara figs;
%     ini_check4artefacts; % Initialisation of variables
%     epochlength = 5;
%     amplitude_threshold_mag = 5000;
%     amplitude_threshold_gra = 10000;
%     out_dir = fullfile(megpath,'ChanCheck');
%     addfiff(sssfile2, strcat(rfxrootdir, 'MF/') ); %not the trans file
%     addtask('maxminamp');   % maximum minus minimum amplitude within epoch
%     addtask('maxmingrad');  % maximum signal gradient (amplitude difference between successive data points in time) with epoch
%     addtask('threshold');   % number of epochs with max-min amplitudes above specified threshold
%     addtask('crosscorr');   % intercorrelation matrix for magnetometers and gradiometers separately
%     addtask('trendamp');    % linear trend in max-min differences across all epochs
%     addtask('trendgrad');   % linear trend in maximum signal gradients across all epochs
%     check4artefacts;    % That's where it all happens... lean back and enjoy!
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% % % ********plot*********
% % % 1) Get the digitised head points from the fif file (in meters in head coordinates)
% %
% % [co,ki, nu] = hpipoints('vstm2_raw.fif');
% %
% % %2) Get the sphere details you want to plot, e.g. from the maxfilter output, or from fit_sphere_to_points. Here I'll just use an arbitrary guess (again, in meters in head coordinates):
% %
% % o=[0 0 0.045]; %(x y z coordinates of origin)
% % r=0.090 %(radius of sphere)
% %
% % %3) Plot it. This is just an example, it can be prettified in whatever way you want:
% %
% % f=spm_figure;
% % set(f,'renderer','opengl')
% % [x y z]=sphere(30);
% % for a=1:3
% %     subplot(2,2,a)
% %     plot3(co(1,:),co(2,:),co(3,:),'r.');
% %     hold on
% %     s=surf(x*r+o(1),y*r+o(2),z*r+o(3),'facecolor','k','facealpha',0.1,'edgecolor','none','facelighting','phong','specularexponent',0.4,'specularstrength',2);
% %     plot3(o(1),o(2),o(3),'k*','markersize',10)
% %     axis tight image
% %     view(90*double(a==2),90*double(a==1))
% %     L=camlight('headlight');
% % end
% 
% 
% 
% 
