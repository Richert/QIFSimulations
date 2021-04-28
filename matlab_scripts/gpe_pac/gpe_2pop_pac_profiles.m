%%
clc;
clear all
close all

%% load input signal and define parameters

% data loading
all_data = load('ss_data.mat');
fnames = fieldnames(all_data);

% data characteristics
N = numel(fnames);
sr = 1000;
cutoff = ceil(1*sr);
cutoff2 = ceil(20*sr);
data_idx = 2;
n_bins = 32;

% filter characteristics
phase_band = [2, 32];
amp_band = [50, 250];
n_low = 2;
n_high = 10;

%% test plotting

data_tmp = all_data.circuit_8.data(data_idx, cutoff:cutoff2);
figure();
plot(data_tmp');

%% calculate PAC and PPC profiles for each condition

stim_freqs = zeros(1,N);
stim_amps = zeros(size(stim_freqs));
PAC_profiles = cell(N, 1);
PPC_profiles = cell(N,1);
PAC_phase_freqs = cell(N, 1);
PAC_amp_freqs = cell(N, 1);

for i=1:N
    
    % extract data and stimulation condition
    fn = fnames(i);
    condition = getfield(all_data, fn{1});
    stim_freqs(i) = condition.omega;
    stim_amps(i) = condition.alpha;
    data = condition.data(data_idx, :);
    
    % calculate PAC and PPC
    [pac, ppc, phase_freqs, amp_freqs] = crossfreqcoupling(data,sr,phase_band,amp_band,n_bins,n_low,n_high);
    PAC_profiles{i} = pac;
    PPC_profiles{i} = ppc;
    PAC_phase_freqs{i} = phase_freqs;
    PAC_amp_freqs{i} = amp_freqs;
    
end

%% alternatively load existing PAC calculations

% load('PAC_oscillatory.mat')
% phase_freqs = PAC_phase_freqs{1};
% amp_freqs = PAC_amp_freqs{1};

%% plot PAC and PPC profile for a number of conditions

condition_indices = [523, 536];
for i=condition_indices

    display(sprintf('\\omega = %.1f, \\alpha = %.1f', stim_freqs(i), stim_amps(i)))
    
    % PAC
    figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', [10 10 900 300]);
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PAC_profiles{i}');
    axis xy;
    colorbar;
    title('PAC', 'FontSize', 16);
    xlabel('phase frequency', 'FontSize', 16)
    ylabel('amplitude frequency', 'FontSize', 16)
    saveas(gcf, sprintf('PAC_ss_%u', i), 'svg')
    
    % PPC
    figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold', 'Position', [10 10 900 300]);
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PPC_profiles{i}');
    axis xy;
    colorbar;
    title('PPC', 'FontSize', 16);
    xlabel('phase frequency', 'FontSize', 16)
    ylabel('amplitude frequency', 'FontSize', 16)
    saveas(gcf, sprintf('PPC_ss_%u', i), 'svg')
    
end


%% save final data to mat files

name = 'PAC_PPC_profiles_ss.mat';
save(name, 'PAC_profiles', 'PAC_phase_freqs', 'PAC_amp_freqs', 'stim_freqs', 'stim_amps', 'PPC_profiles')

