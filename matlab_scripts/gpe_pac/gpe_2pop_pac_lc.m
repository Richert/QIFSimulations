%%
clc;
clear all
close all

%% load input signal and define parameters

% data loading
all_data = load('lc_data.mat');
fnames = fieldnames(all_data);

% data characteristics
N = numel(fnames);
sr = 1000;
cutoff = ceil(1*sr);
cutoff2 = ceil(20*sr);
data_idx = 2;
n_bins = 16;

% filter characteristics
phase_band = [2, 20];
amp_band = [50, 250];
n_window = 1024;
n_overlap = 512;
eval_freqs = phase_band(1):1:amp_band(2);

%% calculate PAC profiles for each condition

stim_freqs = zeros(1,N);
stim_amps = zeros(size(stim_freqs));
PAC_vals = zeros(size(stim_freqs));
PPC_vals = zeros(size(stim_freqs));
PAC_phase_freqs = zeros(size(stim_freqs));
PAC_amp_freqs = zeros(size(stim_freqs));
PSD_vals = zeros(N,length(eval_freqs));

for i=1:N
    fn = fnames(i);
    condition = getfield(all_data, fn{1});
    stim_freqs(i) = condition.omega;
    stim_amps(i) = condition.alpha;
    data = condition.data(data_idx, :);
    [pac, ppc, fphase, famp, ps] = get_cfc_lohi(data,sr,phase_band,amp_band,n_bins,n_window,n_overlap,eval_freqs);
    PAC_vals(i) = pac;
    PPC_vals(i) = ppc;
    PAC_phase_freqs(i) = fphase;
    PAC_amp_freqs(i) = famp;
    PSD_vals(i,:) = ps;
end

%% alternatively load existing PAC calculations

% load('PAC_oscillatory.mat')
% phase_freqs = PAC_phase_freqs{1};
% amp_freqs = PAC_amp_freqs{1};

%% create PAC/PPC matrices over phase and amplitude frequencies

stim_freqs_unique = unique(stim_freqs);
stim_amps_unique = unique(stim_amps);
n_cols = length(stim_freqs_unique);
n_rows = length(stim_amps_unique);
target_freq = stim_freqs_unique(18);

PAC = zeros(n_rows, n_cols);
PAC_fphase = zeros(size(PAC));
PAC_famp = zeros(size(PAC));
PPC = zeros(size(PAC));
PSDs = zeros(n_rows, length(eval_freqs));

for i=1:N
    
    % get the stimulation frequency and amplitude for condition i
    omega = stim_freqs(i);
    alpha = stim_amps(i);
    
    % get the row and column for this combination of stim freq and amp 
    col = find(stim_freqs_unique == omega);
    row = find(stim_amps_unique == alpha);
    
    % save maximum PAC value, and the frequencies at maximum PAC
    PAC(row,col) = PAC_vals(i);
    PPC(row,col) = PPC_vals(i);
    fphase = PAC_phase_freqs(i);
    famp = PAC_amp_freqs(i);
    PAC_fphase(row,col) = fphase;
    PAC_famp(row,col) = famp;
    
    % save PSD profile for a target frequency
    if omega == target_freq
        PSDs(row,:) = PSD_vals(i, :);
    end
    
end

%% Results plotting

% maximum PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC);
axis xy;
colorbar;
saveas(gcf, 'PAC_lohi_lc', 'svg')

% phase frequency at PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_fphase);
axis xy;
colorbar;
saveas(gcf, 'PAC_fphase_lc', 'svg')

% amplitude frequency at PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_famp);
axis xy;
colorbar;
saveas(gcf, 'PAC_famp_lc', 'svg')

% coherence between low and high frequency
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PPC);
axis xy;
colorbar;
saveas(gcf, 'PPC_lohi_lc', 'svg')

% plot PSDs for all amplitudes of a given driving frequency
figure;
imagesc(eval_freqs,stim_amps_unique,log(PSDs))
axis xy;
colorbar;
saveas(gcf, 'PSD_lc', 'svg')

%% save final data to mat files

name = 'CFP_lc.mat';
save(name, 'PAC_vals', 'PAC_phase_freqs', 'PAC_amp_freqs', 'stim_freqs', 'stim_amps', ...
     'stim_freqs_unique', 'stim_amps_unique', 'PAC', 'PAC_fphase', 'PAC_famp', ...
     'PPC')
