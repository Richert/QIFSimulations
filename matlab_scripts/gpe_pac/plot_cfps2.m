clear all;
close all

% load data from file
%%%%%%%%%%%%%%%%%%%%%

path = '/home/rgast/ownCloud/data/gpe_2pop_forced';
condition = 'lc';

fname = sprintf('%s/PAC_PPC_plotting_data_%s.mat', path, condition);
load(fname)

% plot data
%%%%%%%%%%%

plot_size = [10 10 500 400];
font_size1 = 14;
font_size2 = 16;
clim = [0.0, max(PAC_mean(:))];

% mean PAC value
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_mean,clim);
title('mean PAC', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
axis xy;
colorbar;
saveas(gcf, sprintf('PAC_mean_%s', condition), 'svg')

% PAC PPC correlation
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_PPC_corr,[-max(PAC_PPC_corr(:)), max(PAC_PPC_corr(:))]);
title('correlation(PAC,PPC)', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
axis xy;
colorbar;
saveas(gcf, sprintf('PAC_PPC_corr_%s', condition), 'svg')

% mean of PAC*PPC
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_mean_pl,clim);
title('mean PAC * PPC', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
axis xy;
colorbar;
saveas(gcf, sprintf('PAC_mean_pl_%s', condition), 'svg')

% mean of PAC*(1-PPC)
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_mean_npl,clim);
title('mean PAC * (1-PPC)', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
axis xy;
colorbar;
saveas(gcf, sprintf('PAC_mean_npl_%s', condition), 'svg')
