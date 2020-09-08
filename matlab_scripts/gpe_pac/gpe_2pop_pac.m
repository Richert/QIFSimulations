%%
clc;
clear all
close all

%% load input signal and define parameters

% data loading
all_data = load('bs_data.mat');
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

%% calculate PAC profiles for each condition

stim_freqs = zeros(1,N);
stim_amps = zeros(size(stim_freqs));
PAC_profiles = cell(N, 1);
PAC_phase_freqs = cell(N, 1);
PAC_amp_freqs = cell(N, 1);

for i=1:N
    fn = fnames(i);
    condition = getfield(all_data, fn{1});
    stim_freqs(i) = condition.omega;
    stim_amps(i) = condition.alpha;
    data = condition.data(data_idx, :);
    [pac, phase_freqs, amp_freqs] = get_pac_profile(data,sr,phase_band,amp_band,n_bins,n_low,n_high);
    PAC_profiles{i} = pac;
    PAC_phase_freqs{i} = phase_freqs;
    PAC_amp_freqs{i} = amp_freqs;
end

%% alternatively load existing PAC calculations

% load('PAC_oscillatory.mat')
% phase_freqs = PAC_phase_freqs{1};
% amp_freqs = PAC_amp_freqs{1};

%% plot PAC profile for a number of conditions

condition_indices = [10, 20, 100, 200];
for i=condition_indices
    figure(i);
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PAC_profiles{i}');
    axis xy;
    colorbar;
end

%% extract maximum/mean PAC value and its frequency pair for each condition

stim_freqs_unique = unique(stim_freqs);
stim_amps_unique = unique(stim_amps);
n_cols = length(stim_freqs_unique);
n_rows = length(stim_amps_unique);

PAC_max = zeros(n_rows, n_cols);
PAC_mean = zeros(size(PAC_max));
PAC_fphase = zeros(size(PAC_max));
PAC_famp = zeros(size(PAC_max));
PAA_osc = zeros(size(PAC_max));
PAA_env = zeros(size(PAC_max));
n_pac_cols = length(phase_freqs);
n_pac_rows = length(amp_freqs);

for i=1:N
    
    % get the stimulation frequency and amplitude for condition i
    omega = stim_freqs(i);
    alpha = stim_amps(i);
    
    % get the row and column for this combination of stim freq and amp 
    col = find(stim_freqs_unique == omega);
    row = find(stim_amps_unique == alpha);
    
    % save maximum PAC value, and the frequencies at maximum PAC
    PAC_max(row,col) = max(PAC_profiles{i},[],'all');
    PAC_mean(row,col) = mean(PAC_profiles{i},'all');
    [max_row, max_col] = find(PAC_profiles{i}' == PAC_max(row,col));
    fphase = PAC_phase_freqs{i};
    famp = PAC_amp_freqs{i};
    fphase_atmax = fphase(max_col);
    famp_atmax = famp(max_row);
    PAC_fphase(row,col) = fphase_atmax;
    PAC_famp(row,col) = famp(max_row);
    
    % calculate maximum of averaged waveform/envelope at driving frequency
    fn = fnames(i);
    condition = getfield(all_data, fn{1});
    data = condition.data(data_idx, :);
    driver = condition.data(1, :);
    [PLA_po_all, PLA_env_all, PLA_po_max, PLA_env_max] = PhaseLockAmp(data',driver',[stim_freqs(i)-0.5, stim_freqs(i)+0.5],sr,1);
    PAA_osc(row,col) = PLA_po_max;
    PAA_env(row,col) = PLA_env_max;
    
end

%% Results plotting

% maximum PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_max);
axis xy;
colorbar;
saveas(gcf, 'PAC_max_bs', 'svg')

% mean PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_mean);
axis xy;
colorbar;
saveas(gcf, 'PAC_mean_bs', 'svg')


% phase frequency at maximum PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_fphase);
axis xy;
colorbar;
saveas(gcf, 'PAC_fphase_bs', 'svg')

% amplitude frequency at maximum PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_famp);
axis xy;
colorbar;
saveas(gcf, 'PAC_famp_bs', 'svg')

% phase-averaged waveform of observed signal (based on driver phase)
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAA_osc);
axis xy;
colorbar;
saveas(gcf, 'PAA_osc_bs', 'svg')

% phase-averaged envelope of observed signal (based on driver phase)
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAA_osc ./ PAA_env);
axis xy;
colorbar;
saveas(gcf, 'PAA_env_bs', 'svg')

%% save final data to mat files

name = 'PAC_bistable.mat';
save(name, 'PAC_profiles', 'PAC_phase_freqs', 'PAC_amp_freqs', 'stim_freqs', 'stim_amps', ...
     'stim_freqs_unique', 'stim_amps_unique', 'PAC_max', 'PAC_fphase', 'PAC_famp', ...
     'PAA_osc', 'PAA_env')
