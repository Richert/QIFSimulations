function PLV = get_plv(data,srate,low_band,high_band)
%GET_PPC_PROFILE calculate phase-locking value between signal components
%at a low and a high frequency.

low_freq = mean(low_band);
high_freq = mean(high_band);
low_diff = low_band(2)-low_band(1);
high_diff = high_band(2)-high_band(1);
cutoff = ceil(2*srate/low_band(1));

% filter out the low-frequency component
low_filt_order = round(3*(srate/low_freq));
low_filterweights = fir1(low_filt_order,[low_freq-low_diff low_freq+low_diff]/(srate/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
lowfilt = filtfilt(low_filterweights,1,data);
data_lf = lowfilt(cutoff:end-cutoff);

% filter out the low-frequency component
high_filt_order = round(3*(srate/high_freq));
high_filterweights = fir1(high_filt_order,[high_freq-high_diff high_freq+high_diff]/(srate/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
highfilt = filtfilt(high_filterweights,1,data);
data_hf = highfilt(cutoff:end-cutoff);

% get instantaneous phase from analytic signal of filtered components
phase_lf = angle(hilbert(data_lf));
phase_hf = angle(hilbert(data_hf));

% calculate the phase-locking value
[m,n] = rat(low_freq/high_freq);
delta_phi = mod(m*phase_hf-n*phase_lf,2*pi);
PhaseDifference = exp(1i*delta_phi);
PLV = abs(mean(PhaseDifference));

end

