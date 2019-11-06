function movement_figure(subjname,k)

close all
% uses check_movecomp to read and display maxfilter -movecomp output
% Files with maxfilter outputs, e.g. for different subjects and conditions


movecompfiles = {sprintf('%s_%d_ssslogfile.log',subjname,k)};
nr_files = length(movecompfiles);

for ff = 1:nr_files,
    fprintf(1, 'Processing %s\n', movecompfiles{ff});
    [mv_fig, linet, linee, lineg, linev, liner, lined] = check_movecomp(movecompfiles{ff}); % read info from log-file
    set(gcf,'color','w')
    [a,b] = fileparts( movecompfiles{ff} );
%    tittxt = [a(end-17:end) '/' b]; % just one way to keep the figure title short
    figure( mv_fig );
    % check_movecomp only creates the figure, but without the title, so you can do this separately
    ht = title(movecompfiles); set(ht, 'interpreter', 'none'); % add title to figure
end;

[a,b,c]=fileparts(movecompfiles{1});
saveas(gcf,b)