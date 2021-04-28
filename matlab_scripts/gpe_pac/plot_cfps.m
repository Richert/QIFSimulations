clear all;
close all

% load data from file
%%%%%%%%%%%%%%%%%%%%%

path = '/home/rgast/ownCloud/data/gpe_2pop_forced';
condition = 'bs';
condition_indices = [];

fname = sprintf('%s/PAC_PPC_profiles_%s.mat', path, condition);
load(fname)

stim_freqs_unique = sort(unique(stim_freqs));
stim_amps_unique = sort(unique(stim_amps));

% plot data
%%%%%%%%%%%

plot_size = [10 10 500 400];
font_size1 = 14;
font_size2 = 16;
clim = [0.0, 0.16];

% maximum PAC value
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_max,clim);
title('maximum PAC', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
axis xy;
colorbar;
saveas(gcf, sprintf('PAC_max_%s', condition), 'svg')

% maximum PAC times PPC at maximum PAC
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_max.*PPC,clim);
title('PAC * PPC', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
axis xy;
colorbar;
saveas(gcf, sprintf('PAC_and_PPC_%s', condition), 'svg')

% maximum PAC times 1-PPC at maximum PAC
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_max.*(1-PPC),clim);
title('PAC * (1-PPC)', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
axis xy;
colorbar;
saveas(gcf, sprintf('PAC_and_not_PPC_%s', condition), 'svg')

% phase frequency at PAC value
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_fphase);
axis xy;
title('Phase Frequency', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
colorbar;
saveas(gcf, sprintf('PAC_fphase_%s', condition), 'svg')

% amplitude frequency at PAC value
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_famp);
axis xy;
title('Amplitude Frequency', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
colorbar;
saveas(gcf, sprintf('PAC_famp_%s', condition), 'svg')

% coherence between low and high frequency
figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PPC);
axis xy;
title('PPC at maximum PAC', 'FontSize', font_size2);
xlabel('\omega', 'FontSize', font_size2)
ylabel('\alpha', 'FontSize', font_size2)
colorbar;
saveas(gcf, sprintf('PPC_atmax_%s', condition), 'svg')

% specific PAC/PPC profiles
for i=condition_indices
    
    msg = sprintf('\\omega = %.1f, \\alpha = %.1f', stim_freqs(i), stim_amps(i));
    display(msg)
    
    % PAC
    figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PAC_profiles{i}');
    axis xy;
    colorbar;
    title('PAC', 'FontSize', font_size2);
    xlabel('phase frequency', 'FontSize', font_size2)
    ylabel('amplitude frequency', 'FontSize', font_size2)
    saveas(gcf, sprintf('PAC_%s_%u', condition, i), 'svg')
    
    % PPC
    figure('DefaultAxesFontSize', font_size1, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PPC_profiles{i}');
    axis xy;
    colorbar;
    title('PPC', 'FontSize', font_size2);
    xlabel('phase frequency', 'FontSize', font_size2)
    ylabel('amplitude frequency', 'FontSize', font_size2)
    saveas(gcf, sprintf('PPC_%s_%u', condition, i), 'svg')
    
end
