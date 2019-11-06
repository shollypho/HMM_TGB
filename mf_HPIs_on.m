% Record Translation and Rotation movements to assess the amount of
% movement throughout the resting block


all_data = data_setup;

% path to all maxfiltered files
mf_path = '/imaging/hp02/TGB/rest_closed/MF';
cd(mf_path)
cd Patients
result = [];

for ss = 1:length(all_data) % subjects
    
    mf_path = '/imaging/hp02/TGB/rest_closed/MF';
    cd(mf_path)
    if strcmp(all_data{ss,2},'Con')
        cd Controls
    else
        cd Patients
    end
    
    cd(all_data{ss,3})
    
    
    % check if trans file exists
    if exist(sprintf('movecomp_%s_%d_ssslogfile.log',all_data{ss,3},all_data{ss,5}), 'file')
        
        [linet, lineg, linev, liner, lined] = record_movecomp(sprintf('movecomp_%s_%d_ssslogfile.log',all_data{ss,3},all_data{ss,5}));
        hpi_start_250hz(ss) = (length(lineg)-sum(lineg>0.99))*1000/4;
        
    end

end

