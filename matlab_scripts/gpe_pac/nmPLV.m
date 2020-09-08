%% phase-phase coupling
clear;
clc;

%% load data
all_data = load('PAC_bistable');
data = all_data.circuit_51.data(2,:);
data = data(10000:end);
%data = data+randn(size(data)).*600.0; %add noise
fs = 10000;
ts = (0:1:(length(data)-1))/fs;
% plot raw data
figure;
plot(ts(1:10000),data(1:10000));

%% fft calculation
% Fourier = fft(data);
% P2 = abs(Fourier/length(data));
% P1 = P2(1:length(data)/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(length(data)/2))/length(data);
% figure;
% plot(f,P1)
% xlim([2,10]);

%% ------------------------- single pair---------------------%%
%% power spectrum calculation
window = 2*fs;
overlap = 0.5;
nfft = window;
Pxx = zeros(1,nfft/2+1);
Fp = zeros(1,nfft/2+1);
[Pxx,Fp] = pwelch(data,window,[],nfft,fs);

figure;
plot(Fp, Pxx);
xlim([1,200]);

%% find center low frequency and high frequency
ind_low = find(Fp>=1 & Fp<=30);
[peak_low,loc_low] = max(Pxx(ind_low));
freq_low = Fp(ind_low(loc_low));

ind_high = find(Fp>freq_low*2);
[peak_high,loc_high] = max(Pxx(ind_high));
freq_high = Fp(ind_high(loc_high));

%% bandpass filtering
lowband = [freq_low-1,freq_low+1];
low_filt_order    = round(3*(fs/lowband(1)));
low_filterweights = fir1(low_filt_order,[mean(lowband)-(lowband(2)-lowband(1)) mean(lowband)+(lowband(2)-lowband(1))]/(fs/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path

highband = [freq_high-freq_low/2,freq_high+freq_low/2];
high_filt_order    = round(3*(fs/highband(1)));
high_filterweights = fir1(high_filt_order,[mean(highband)-(highband(2)-highband(1)) mean(highband)+(highband(2)-highband(1))]/(fs/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path

lowfilt = filtfilt(low_filterweights,1,data);
highfilt = filtfilt(high_filterweights,1,data);

% -- compute low phase angles and high phase angles
% -- first and last 10 cycls of lowest frequency are removed to account for edge artifacts
margin = 10*ceil(fs*1/(freq_low-1));
lowphase = angle(hilbert(lowfilt((1+margin):end-margin)));
highphase = angle(hilbert(highfilt((1+margin):end-margin)));

ts = ts((1+margin):end-margin);
figure;
plot(ts,lowphase);
hold on;
plot(ts,highphase);
xlim([1,2])

%% n:m phase locking value
% g = gcd(freq_low,freq_high);
% p = freq_low/g;
% q = freq_high/g;
[N,D] = rat(freq_low/freq_high);
m = N;
n = D;
delta_phi = mod(m*highphase-n*lowphase,2*pi);
PhaseDifference = exp(1i*delta_phi);
PLV = abs(mean(PhaseDifference));

%% -------------------------multiple pair---------------------- %%
% frequency range
pha_locut = 3 ;
pha_hicut = 31;
num_phabins = (pha_hicut-pha_locut)/2;

amp_locut = 3;
amp_hicut = 101;
num_ampbins = (amp_hicut-amp_locut)/2;

fphase = linspace(3,31,num_phabins+1);
famp = linspace(3,101,num_ampbins+1);

fphaseCenters = fphase(2:end)-diff(fphase)/2;
fampCenters = famp(2:end)-diff(famp)/2;

PLV = zeros(num_phabins,num_ampbins);

lowphase = cell(num_phabins,1);
highphase = cell(num_phabins,1);

%% filt phase from 4 to 50Hz (2Hz bandwidth)

for idx = 1:num_phabins
    lowfilt = [];
    Hdata_low = [];
    %lowband = [fphase(idx),fphase(idx+1)];
    low_filt_order    = round(3*(fs/fphase(idx)));
    low_filterweights = fir1(low_filt_order,[fphase(idx) fphase(idx+1)]/(fs/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
    lowfilt = filtfilt(low_filterweights,1,data);
    margin = 10*ceil(fs*1/fphase(idx));
    Hdata_low = hilbert(lowfilt((1+margin):end-margin));
    lowphase{idx} = angle(Hdata_low);
end

%% filt amplitude from 4 to 150Hz (4Hz bandwidth)
for indlow = 1:num_phabins
    Hdata_high{indlow} = [];
    margin = 10*ceil(fs*1/fphase(indlow));
    start_amp = min(find(fampCenters > fphaseCenters(indlow)*3/2));
    for idx = start_amp:num_ampbins
        highband = [];
        highband = [(fampCenters(idx)-fphaseCenters(indlow)/2),(fampCenters(idx)+fphaseCenters(indlow)/2)];
        high_filt_order    = round(3*(fs/highband(1)));
        high_filterweights = fir1(high_filt_order,[highband(1) highband(2)]/(fs/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
        highfilt = filtfilt(high_filterweights,1,data);
        Hdata_high{indlow}(idx,:) = hilbert(highfilt((1+margin):end-margin));
    end
    highphase{indlow} = angle(Hdata_high{indlow});
end

%% n:m phase locking value

for indlow = 1:num_phabins
    freq_low = fphaseCenters(indlow);
    start_amp = min(find(fampCenters > fphaseCenters(indlow)*3/2));
    for indhigh = start_amp:num_ampbins
        freq_high = fampCenters(indhigh);
        [N,D] = rat(freq_low/freq_high);
        m = N;
        n = D;
        delta_phi = mod(m*highphase{indlow}(indhigh,:)-n*lowphase{indlow},2*pi);
        PhaseDifference = exp(1i*delta_phi);
        PLV(indlow,indhigh) = abs(mean(PhaseDifference));
    end
end

figure;
imagesc(fphaseCenters,fampCenters,PLV');
axis xy;
colorbar;


%% ------- loop in different alpha setting (or Omega setting) ------- %%
clear;
clc;
fs = 10000;
%% load data
dataset = load('/data/pt_02038/Simulating_data/PAC_bistable');
names = fieldnames(dataset);

subStr = 'circuit_6\d'; %subStr = 'circuit_\d6' for Omega selection
filteredStruct = names(find(~cellfun(@isempty,regexp(names,subStr))));
n_alpha = length(filteredStruct); 
amp_locut = 3;
amp_hicut = 101;
num_ampbins = (amp_hicut-amp_locut)/2;
famp = linspace(3,101,num_ampbins+1);
fampCenters = famp(2:end)-diff(famp)/2;

PLV = zeros(n_alpha,num_ampbins);
alpha = zeros(n_alpha,1);
for i = 1:n_alpha
    data = dataset.(filteredStruct{i+1}).data(2,:);
    data = data(10000:end);
    %data = data+randn(size(data)).*600.0; %add noise
    ts = (0:1:(length(data)-1))/fs;
    freq_low = dataset.(filteredStruct{i+1}).omega;
    alpha(i) = dataset.(filteredStruct{i+1}).alpha;
    %% bandpass filtering
    lowband = [freq_low-1,freq_low+1];
    low_filt_order    = round(3*(fs/lowband(1)));
    low_filterweights = fir1(low_filt_order,[lowband(1) lowband(2)]/(fs/2));
    lowfilt = filtfilt(low_filterweights,1,data);
    margin = 10*ceil(fs*1/(freq_low-1));
    lowphase = angle(hilbert(lowfilt((1+margin):end-margin)));
    
    start_amp = min(find(fampCenters > freq_low*3/2));
    for idx = start_amp:num_ampbins
        highband = [];
        highphase = [];
        highband = [(fampCenters(idx)-freq_low/2),(fampCenters(idx)+freq_low/2)];
        high_filt_order    = round(3*(fs/highband(1)));
        high_filterweights = fir1(high_filt_order,[highband(1) highband(2)]/(fs/2)); % -- NOTE: if eeglab is added to path, remove eeglab/plugins/biosig from path
        highfilt = filtfilt(high_filterweights,1,data);
        highphase = angle(hilbert(highfilt((1+margin):end-margin)));
        %% calculating the PLV
        freq_high = fampCenters(idx);
        [N,D] = rat(freq_low/freq_high);
        m = N;
        n = D;
        delta_phi = mod(m*highphase-n*lowphase,2*pi);
        PhaseDifference = exp(1i*delta_phi);
        PLV(i,idx) = abs(mean(PhaseDifference));
    end    
end

f=figure;
f.Position= [348,175,449,437];
imagesc(alpha,fampCenters,PLV');
axis xy;
colorbar;