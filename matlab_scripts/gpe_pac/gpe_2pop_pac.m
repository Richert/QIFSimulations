%%
clc;
clear all
close all

%% load existing PAC calculations

condition = 'lc';
load(sprintf('/home/rgast/ownCloud/data/gpe_2pop_forced/PAC_PPC_profiles_%s.mat', condition));
phase_freqs = PAC_phase_freqs{1};
amp_freqs = PAC_amp_freqs{1};

%% plot PAC and PPC profile for a number of conditions

condition_indices = [10, 20, 100, 120];
for i=condition_indices
    
    % PAC
    figure();
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PAC_profiles{i}');
    axis xy;
    colorbar;
    title(sprintf('PAC: \\alpha = %.2f, \\omega = %.2f', stim_amps(i), stim_freqs(i)))
    
    figure();
    imagesc(PAC_phase_freqs{i},PAC_amp_freqs{i},PPC_profiles{i}');
    axis xy;
    colorbar;
    title(sprintf('PPC: \\alpha = %.2f, \\omega = %.2f', stim_amps(i), stim_freqs(i)))
end

%% extract maximum/mean PAC value and its frequency pair for each condition

stim_freqs_unique = unique(stim_freqs);
stim_amps_unique = unique(stim_amps);
n_cols = length(stim_freqs_unique);
n_rows = length(stim_amps_unique);

PAC_max = zeros(n_rows, n_cols);
PAC_max_pl = zeros(size(PAC_max));
PAC_max_npl = zeros(size(PAC_max));
PAC_mean = zeros(size(PAC_max));
PAC_mean_pl = zeros(size(PAC_max));
PAC_mean_npl = zeros(size(PAC_max));

PAC_fphase = zeros(size(PAC_max));
PAC_famp = zeros(size(PAC_max));

PPC_atmax = zeros(size(PAC_max));
PPC_max = zeros(size(PAC_max));
PPC_mean = zeros(size(PAC_max));

PAC_PPC_corr = zeros(size(PAC_max));

for i=1:length(PAC_profiles)
    
    % get data for condition i
    pac = PAC_profiles{i}';
    ppc = PPC_profiles{i}';
    omega = stim_freqs(i);
    alpha = stim_amps(i);
    
    % get the row and column for this combination of stim freq and amp 
    col = find(stim_freqs_unique == omega);
    row = find(stim_amps_unique == alpha);
    
    % save PAC aggregates
    PAC_max(row,col) = max(pac,[],'all');
    PAC_max_pl(row,col) = max(pac.*ppc,[],'all');
    PAC_max_npl(row,col) = max(pac.*(1-ppc),[],'all');
    PAC_mean(row,col) = mean(pac,'all');
    PAC_mean_pl(row,col) = mean(pac.*ppc, 'all');
    PAC_mean_npl(row,col) = mean(pac.*(1-ppc), 'all');
    
    % save phase freq, amp freq and ppc at max PAC
    [max_row, max_col] = find(pac == PAC_max(row,col));
    fphase = PAC_phase_freqs{i};
    famp = PAC_amp_freqs{i};
    fphase_atmax = fphase(max_col);
    famp_atmax = famp(max_row);
    PAC_fphase(row,col) = fphase_atmax;
    PAC_famp(row,col) = famp(max_row);
    PPC_atmax(row,col) = ppc(max_row,max_col);
    PPC_max(row,col) = max(ppc,[],'all');
    PPC_mean(row,col) = mean(ppc,'all');
    
    % calculate correlation between PAC and PPC profiles
    pac_z = zscore(pac(:));
    ppc_z = zscore(ppc(:));
    c = corrcoef(pac_z, ppc_z);
    PAC_PPC_corr(row,col) = c(1,2);
    
end

%% Results plotting

% maximum PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_max);
axis xy;
colorbar;
title('max PAC')

% maximum PAC*PPC value 
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_max_pl);
axis xy;
colorbar;
title('max PAC*PPC')

% maximum PAC*(1-PPC) value 
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_max_npl);
axis xy;
colorbar;
title('max PAC*(1-PPC)')

% mean PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_mean);
axis xy;
colorbar;
title('mean PAC')

% mean PAC*PPC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_mean_pl);
axis xy;
colorbar;
title('mean PAC*PPC')

% mean PAC*(1-PPC) value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_mean_npl);
axis xy;
colorbar;
title('mean PAC*(1-PPC)')

% phase frequency at maximum PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_fphase);
axis xy;
colorbar;
title('Frequency of phase at max PAC')

% amplitude frequency at maximum PAC value
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_famp);
axis xy;
colorbar;
title('Frequency of amplitude at max PAC')

% PPC at max PAC
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PPC_atmax);
axis xy;
colorbar;
title('PPC at max PAC')

% max PPC
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PPC_max);
axis xy;
colorbar;
title('max PPC')

% mean PPC
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PPC_max);
axis xy;
colorbar;
title('mean PPC')

% PAC-PPC correlation
figure;
imagesc(stim_freqs_unique,stim_amps_unique,PAC_PPC_corr);
axis xy;
colorbar;
title('corr(PAC,PPC)')

%% save final data to mat files

name = sprintf('PAC_PPC_plotting_data_%s.mat', condition);
save(name, 'PAC_max', 'PAC_max_pl', 'PAC_max_npl', 'PAC_mean', 'PAC_mean_pl', ...
    'PAC_mean_npl', 'PPC_atmax', 'PPC_max', 'PPC_mean', 'PAC_famp', 'PAC_fphase', ...
    'stim_freqs_unique', 'stim_amps_unique', 'PAC_PPC_corr')
