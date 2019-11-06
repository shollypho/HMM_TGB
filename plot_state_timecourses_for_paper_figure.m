addpath(genpath('/imaging/hp02/software_n_scripts/linspecer'));

colr = {'r', 'b', 'k', 'g', 'c', 'm'};

gamma2 = Gamma;
gamma2(gamma2>0.5) = 1;
gamma2(gamma2<=0.5) = 0;

figure('Color','w','position',[20 72 1200 600])
K = 12;
C = linspecer(K);
subplot(2,1,1)
for i= 1:6
    plot(1:1500, (gamma2(1:1500,i)+i+(0.5*i)), 'color', C(i*2,:), 'linewidth', 1); 
    ylim([0 11])
    hold on
end

subplot(2,1,2)
for i= 1:6
    plot(1:1500, (data(i,1001:2500)), 'color', C(i*2,:), 'linewidth', 1); 
    %ylim([0 11])
    hold on
end

% Cov matrix
A = data(1:20,1:500)';
cov_mat = cov(A);
figure; imagesc(cov_mat); colormap(linspecer);
A = rand(15);
figure; imagesc(A);
figure; imagesc(A); colormap(linspecer);