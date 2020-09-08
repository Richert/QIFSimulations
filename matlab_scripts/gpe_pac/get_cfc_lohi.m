function [pac, ppc, freq_lo, freq_hi, powers] = get_cfc_lohi(data,srate,phase_band,amp_band,n_bins,n_window,n_overlap,freqs)
%GET_PAC_PROFILE calculates PAC profile of a 1d timeseries (data) with a
% sampling rate (srate) for combinations of phase frequencies and amplitude 
% frequencies within the provided bands. 

% frequency binning 
%%%%%%%%%%%%%%%%%%%

pha_locut = phase_band(1);
pha_hicut = phase_band(2);
cutoff = 2*srate/pha_locut;
pnts = length(data(:, cutoff:end-cutoff));
bw = 8;

amp_locut = amp_band(1);
amp_hicut = amp_band(2);

% find frequencies around which to filter the signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[powers, frequencies] = pwelch(data, hamming(n_window), n_overlap, freqs, srate);
idx_lo = (frequencies >= pha_locut) .* (frequencies <= pha_hicut);
idx_hi = (frequencies >= amp_locut) .* (frequencies <= amp_hicut);
[~, max_idx_lo] = max(powers(idx_lo == 1));
[~, max_idx_hi] = max(powers(idx_hi == 1));
freq_lo = frequencies(idx_lo == 1);
freq_lo = freq_lo(max_idx_lo);
freq_hi = frequencies(idx_hi == 1);
freq_hi = freq_hi(max_idx_hi);

% filtering of the signal around the low and high frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filtering of the phase-containing frequency
lowband = [freq_lo-ceil(freq_lo/bw),freq_lo+ceil(freq_lo/bw)];
low_filt_order = round(3*(srate/lowband(1)));
low_filterweights = fir1(low_filt_order,[mean(lowband)-(lowband(2)-lowband(1)) mean(lowband)+(lowband(2)-lowband(1))]/(srate/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
lowfilt = filtfilt(low_filterweights,1,data);
data4phase = lowfilt(cutoff:end-cutoff);
getphase = angle(hilbert(data4phase));

% filtering of the amplitude-containing frequencies
highband = [freq_hi-ceil(freq_hi/bw),freq_hi+ceil(freq_hi/bw)];
high_filt_order = round(3*(srate/highband(1)));
high_filterweights = fir1(high_filt_order,[mean(highband)-(highband(2)-highband(1)) mean(highband)+(highband(2)-highband(1))]/(srate/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
highfilt = filtfilt(high_filterweights,1,data);
data4amp = highfilt(cutoff:end-cutoff);
getamp = abs(hilbert(data4amp));

% PAC and PPC calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pac = get_pac(getphase,getamp,n_bins);
% [coh, ~] = mscohere(data4phase, data4amp, hamming(n_window), n_overlap, [lowband(1), freq_lo, lowband(2)], srate);
% ppc = coh(2);
[~,~,PLAosci,PLAamp] = PhaseLockAmp(data4amp', data4phase', phase_band, srate, 2);
ppc = PLAosci/PLAamp; 

end