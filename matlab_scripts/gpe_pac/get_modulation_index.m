%% function for calculating the PAC modulation index
function MI = get_modulation_index(raw_phase,raw_amp,n_bins)

% extract information
n_phabins = size(raw_phase,1);
n_ampbins = size(raw_amp{1},1);
n_timepoints = size(raw_amp{1},2);

% calculate PAC (modulation index) for each pair of low- and high frequency
MI = zeros(n_phabins,n_ampbins);

for i = 1:n_phabins
    for j = 1:n_ampbins
        
        % get low freq phase and high freq amp
        low_phase = raw_phase(i,:);
        high_amp = raw_amp{i}(j,:);
        
        lowphase_bin = ceil( tiedrank( low_phase ) / (n_timepoints / n_bins) ); % -- bin the low phase angles into nbins -- NOTE: tiedrank also exists in eeglab toolbox; when added to path, may cause conflict
        highamp_bin = zeros(1,n_bins);
        for k=1:n_bins
            highamp_bin(k) = squeeze(mean(high_amp(lowphase_bin==k))); % -- compute mean high power in each bin
        end
        highamp_bin = highamp_bin ./ sum(highamp_bin); % -- normalize
        MI(i,j) = 1 + sum(highamp_bin.*log(highamp_bin)) / log(n_bins);
    end
end

end
