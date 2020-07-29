%% scripts for Retreat PAC calculation
%% add function
% addpath('D:\study\RetreatScripts\MI_PAC_adjust.m');
%%
clc;
clear;

%% load input signal
all_data = load('PAC_oscillatory.mat');
data = all_data.circuit_0.data(2, 1:end);
%data = data + randn(size(data)).*0.01;
figure();
plot(data');
nbins = 64;

%% single frequency pair
% frequency range
lowband = [1 2];
highband = [120 160];
srate = 10000;

% bandpass filter
low_filt_order = round(3*(srate/lowband(1)));
low_filterweights = fir1(low_filt_order,[mean(lowband)-(lowband(2)-lowband(1)) mean(lowband)+(lowband(2)-lowband(1))]/(srate/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path

high_filt_order = round(3*(srate/highband(1)));
high_filterweights = fir1(high_filt_order,[mean(highband)-(highband(2)-highband(1)) mean(highband)+(highband(2)-highband(1))]/(srate/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path

lowfilt = filtfilt(low_filterweights,1,data);
highfilt = filtfilt(high_filterweights,1,data);

% -- compute low phase angles and high power
% -- first and last 1000 pnts are removed to account for edge artifacts
cutoff = 2*srate/1;
lowphase = angle(hilbert(lowfilt(cutoff:end-cutoff)));
highpow = abs(hilbert(highfilt(cutoff:end-cutoff)));

ntimepoints = length(highpow);
% compute MI
lowphase_bin = ceil( tiedrank( lowphase ) / (ntimepoints / nbins) ); % -- bin the low phase angles into nbins -- NOTE: tiedrank also exists in eeglab toolbox; when added to path, may cause conflict
highpow_bin = zeros(1,nbins);
for k=1:nbins
    highpow_bin(k) = squeeze(mean(highpow(lowphase_bin==k))); % -- compute mean high power in each bin
end
highpow_bin = highpow_bin ./ sum(highpow_bin); % -- normalize
tmpmi = (log(nbins) + sum(highpow_bin.*log(highpow_bin)) ) ./ log(nbins);

% plot MI histogram 
figure;
bar(1:nbins,highpow_bin,'hist')
set(gca,'ylim',[0 0.2],'FontName','Arial','FontSize',10)
%set(gca,'xlim',[0 19])
text(2,0.18,['MI = ' num2str(tmpmi)]);

%% multi- frequency pair
% frequency range
pha_locut = 1;
pha_hicut = 20;
num_phabins = round((pha_hicut-pha_locut)/1);
cutoff = 2*srate/pha_locut;
pnts = length(data(:, cutoff:end-cutoff));

amp_locut = 20;
amp_hicut = 120;
num_ampbins = round((amp_hicut-amp_locut)/10);

fphase = linspace(pha_locut,pha_hicut,num_phabins+1);
famp = linspace(amp_locut,amp_hicut,num_ampbins+1);

fphaseCenters = fphase(2:end)-diff(fphase)/2;
fampCenters = famp(2:end)-diff(famp)/2;

MI = zeros(num_phabins,num_ampbins);

data4phase= [];
data4amp = [];
Hdata4phase = [];
Hdata4amp = [];

data4phase= [];
data4amp = cell(num_phabins,1);
Hdata4phase = [];
Hdata4amp = cell(num_phabins,1);

%% filt phase from 4 to 50Hz (2Hz bandwidth)

for idx = 1:num_phabins
    lowband = [];
    lowband = [fphase(idx),fphase(idx+1)];
    low_filt_order = round(3*(srate/lowband(1)));
    low_filterweights = fir1(low_filt_order,[mean(lowband)-(lowband(2)-lowband(1)) mean(lowband)+(lowband(2)-lowband(1))]/(srate/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
    lowfilt = filtfilt(low_filterweights,1,data);
    data4phase(idx,:) = lowfilt(:, cutoff:end-cutoff);
    Hdata4phase(idx,:) = hilbert(data4phase(idx,:));
end
getphase = angle(Hdata4phase);

%% filt amplitude from 4 to 150Hz (4Hz bandwidth)

for indpha = 1:num_phabins
    data4amp{indpha} = zeros(num_ampbins,pnts);
    Hdata4amp{indpha} = zeros(num_ampbins,pnts);
    start_amp = min(find(fampCenters > fphaseCenters(indpha)*3/2));
    for idx = start_amp:num_ampbins
        highband = [];
        highband = [(fampCenters(idx)-fphaseCenters(indpha)/2),(fampCenters(idx)+fphaseCenters(indpha)/2)];
        high_filt_order = round(3*(srate/highband(1)));
        high_filterweights = fir1(high_filt_order,[mean(highband)-(highband(2)-highband(1)) mean(highband)+(highband(2)-highband(1))]/(srate/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
        highfilt = filtfilt(high_filterweights,1,data);
        data4amp{indpha}(idx,:) = highfilt(cutoff:end-cutoff);
        Hdata4amp{indpha}(idx,:) = hilbert(data4amp{indpha}(idx,:));
        getamp{indpha} = abs(Hdata4amp{indpha});
    end
end

MI(:,:) = MI_PAC_adjust(getphase,getamp,nbins);

% plot MI colormap
figure;
imagesc(fphaseCenters,fampCenters,MI');
axis xy;
colorbar;