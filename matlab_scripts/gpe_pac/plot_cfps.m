clear all;
close all

% load data from file
%%%%%%%%%%%%%%%%%%%%%

condition = 'lc';
condition_indices = [523, 536];
fname1 = sprintf('PAC_PPC_profiles_%s.mat', condition);
fname2 = sprintf('CFP_%s.mat', condition);

load(fname1)
load(fname2)

% plot data
%%%%%%%%%%%

plot_size = [10 10 600 400];

% maximum PAC value
figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC);
title('PAC', 'FontSize', 16);
xlabel('\omega', 'FontSize', 16)
ylabel('\alpha', 'FontSize', 16)
axis xy;
colorbar;
saveas(gcf, 'PAC_lohi_lc', 'svg')

% phase frequency at PAC value
figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_fphase);
axis xy;
title('Phase Frequency', 'FontSize', 16);
xlabel('\omega', 'FontSize', 16)
ylabel('\alpha', 'FontSize', 16)
colorbar;
saveas(gcf, 'PAC_fphase_lc', 'svg')

% amplitude frequency at PAC value
figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PAC_famp);
axis xy;
title('Amplitude Frequency', 'FontSize', 16);
xlabel('\omega', 'FontSize', 16)
ylabel('\alpha', 'FontSize', 16)
colorbar;
saveas(gcf, 'PAC_famp_lc', 'svg')

% coherence between low and high frequency
figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(stim_freqs_unique,stim_amps_unique,PPC);
axis xy;
title('PPC', 'FontSize', 16);
xlabel('\omega', 'FontSize', 16)
ylabel('\alpha', 'FontSize', 16)
colorbar;
saveas(gcf, 'PPC_lohi_lc', 'svg')

% plot PSDs for all amplitudes of a given driving frequency
figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
imagesc(eval_freqs,stim_amps_unique,log(PSDs))
axis xy;
title('PSD', 'FontSize', 16);
xlabel('frequency', 'FontSize', 16)
ylabel('\alpha', 'FontSize', 16)
colorbar;
saveas(gcf, 'PSD_lc', 'svg')

% specific PAC/PPC profiles
for i=condition_indices
    
    msg = sprintf('\\omega = %.1f, \\alpha = %.1f', stim_freqs(i), stim_amps(i));
    display(msg)
    
    % PAC
    figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PAC_profiles{i}');
    axis xy;
    colorbar;
    title('PAC', 'FontSize', 16);
    xlabel('phase frequency', 'FontSize', 16)
    ylabel('amplitude frequency', 'FontSize', 16)
    saveas(gcf, sprintf('PAC_ss_%u', i), 'svg')
    
    % PPC
    figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PPC_profiles{i}');
    axis xy;
    colorbar;
    title('PPC', 'FontSize', 16);
    xlabel('phase frequency', 'FontSize', 16)
    ylabel('amplitude frequency', 'FontSize', 16)
    saveas(gcf, sprintf('PPC_ss_%u', i), 'svg')
    
end
