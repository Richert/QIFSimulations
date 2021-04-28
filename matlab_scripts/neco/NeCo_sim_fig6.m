close all; clear all;

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 1e-3;
dt2 = 1e-2;
N2 = 5e1;
T_sim = 1020.0;
T = round(T_sim/dt);
T_init = 20.0;
T_rec = round((T_sim - T_init)/dt);
T_store = T_init/dt;
v_mac_init = -2;

% QIF parameters
Delta = 2.0;
eta = -5.5;
J = 15.*sqrt(Delta);
v_th = 1e2;

% synpatic efficacy adaptation
tau_a = 10.0;
alpha = 0.05;

% sweep parameters
Ns = [1000, 2000, 4000, 8000];
Ps = [0.01, 0.03, 0.1, 0.3, 1.0];
freqs = zeros(length(Ns), length(Ps));
amps = zeros(size(freqs));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% macroscopic simulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% macroscopic state variables
r_mac_1 = 0.;
v_mac_1 = v_mac_init;
e_mac_1 = 0.;
a_mac_1 = 0.;

% macroscopic observation variables
r_mac_rec = zeros(1,T);
v_mac_rec = zeros(size(r_mac_rec));
e_mac_rec = zeros(size(r_mac_rec));
r_mac_rec_av = zeros(1,T_rec.*dt2);
v_mac_rec_av = zeros(size(r_mac_rec_av));
e_mac_rec_av = zeros(size(r_mac_rec_av));

% input definition
mu = 0.0; % -> for the low activity regime (eta = -6.0)
stim = T_init; % -> for the low activity regime (eta = -6.0)
dur = 200.0;
It = dt.*[1:T];
I = 0.*((It>0)&(It<T_sim)) + mu.*((It>stim)&(It<stim+dur));

for t = 1:T

    % calculate input
    IJ_mac = J.*r_mac_1.*(1 - e_mac_1);
    
    % macroscopic state variable evolution
    r_mac_2 = r_mac_1 + dt .* (Delta./pi + 2 .* v_mac_1.*r_mac_1);
    v_mac_2 = v_mac_1 + dt .* (v_mac_1.^2 + eta + I(t) + IJ_mac - (pi.*r_mac_1).^2);   
    e_mac_2 = e_mac_1 + dt .* a_mac_1;
    a_mac_2 = a_mac_1 + dt .* (alpha.*r_mac_1./tau_a - 2.*a_mac_1./tau_a - e_mac_1./tau_a.^2);
    
    % macroscopic state variable updates
    v_mac_1 = v_mac_2;
    r_mac_1 = r_mac_2;
    a_mac_1 = a_mac_2;
    e_mac_1 = e_mac_2;
    
    % macroscopic state variable recordings
    r_mac_rec(t) = r_mac_2;
    v_mac_rec(t) = v_mac_2;
    e_mac_rec(t) = e_mac_2;
    
    % coarse grained observations
    if t > T_store

        t_tmp = t-T_store;

        if (mod(t,1/dt2)==0)

            % macroscopic state variables
            v_mac_rec_av(t_tmp.*dt2) = mean(v_mac_rec(t-1/dt2+1:t));
            r_mac_rec_av(t_tmp.*dt2) = mean(r_mac_rec(t-1/dt2+1:t));
            e_mac_rec_av(t_tmp.*dt2) = mean(e_mac_rec(t-1/dt2+1:t));

        end
    end
                
end

% macroscopic burst detection
peaks_mac = findpeaks(r_mac_rec_av,'MinPeakProminence',0.5,'MinPeakDistance',100);
freq_mac = length(peaks_mac)/5;
amp_mac = max(r_mac_rec_av(1000:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% microscopic simulations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;
for N=Ns
    j = 1;
    for p=Ps

        % initializations
        %%%%%%%%%%%%%%%%%

        eta_0 = eta+Delta.*tan((pi/2).*(2.*[1:N]-N-1)./(N+1));
        if p < 1.0
            C = rand(N, N);
            C_sorted = sort(C(:), 'ascend');
            C = (C <= C_sorted(N^2*p)) * 1.0;
        end

        % microscopic state variables
        v_1 = v_mac_init.*ones(size(eta_0));
        e_1 = zeros(size(v_1));
        a_1 = zeros(size(v_1));

        % spiking mechanism variables
        wait = zeros(1,N);
        spike_n = zeros(1,N);
        spike_t = wait;

        % microscopic observation variables
        r_rec = zeros(1,T);
        v_rec = zeros(size(r_rec));
        e_rec = zeros(size(r_rec));
        r_rec_av = zeros(1,T_rec.*dt2);
        v_rec_av = zeros(size(r_rec_av));
        e_rec_av = zeros(size(r_rec_av));

        % simulation
        %%%%%%%%%%%%

        for t = 1:T

            % calculate current spikes
            [fr, spike_t, spike_n, mask, wait] = get_network_input(spike_t, wait, dt);
            if p < 1
                s = fr*C./(N*p);
            else
                s = mean(fr);
            end
            I_mic = eta_0 + I(t) + J.*s.*(1 - min(1, e_1));

            % microscopic state variable evolution
            v_2 = v_1 + dt .* mask .* (v_1.^2 + I_mic);
            e_2 = e_1 + dt .* a_1;
            a_2 = a_1 + dt .* (alpha.*s./tau_a - 2.*a_1./tau_a - e_1./tau_a.^2);

            % microscopic spiking mechanism
            [v_2, spike_t, wait] = spiking_mechanism(v_2, v_th, dt, wait, spike_t, 1, I_mic);

            % microscopic state variable updates
            v_1 = v_2;
            a_1 = a_2;
            e_1 = e_2;

            % microscopic state variable recordings
            r_rec(t) = mean(spike_n)./dt;
            v_rec(t) = mean(v_2(wait==0));
            e_rec(t) = mean(e_2(wait==0));

            % coarse grained observations
            if t > T_store

                t_tmp = t-T_store;

                if (mod(t,1/dt2)==0)

                    % microscopic state variables
                    v_rec_av(t_tmp.*dt2) = mean(v_rec(t-1/dt2+1:t));
                    r_rec_av(t_tmp.*dt2) = mean(r_rec(t-1/dt2+1:t));
                    e_rec_av(t_tmp.*dt2) = mean(e_rec(t-1/dt2+1:t));

                end
            end

        end

        % evaluate match between macroscopic and microscopic description
        peaks = findpeaks(r_rec_av,'MinPeakProminence',0.5,'MinPeakDistance',100);
        freqs(i,j) = length(peaks)/5;
        amps(i,j) = max(r_rec_av(1000:end));

        j = j+1;

    end
    i = i+1;
end

save('neco_fig6_data.mat', 'freqs', 'amps', 'freq_mac', 'amp_mac')

figure()        
imagesc(freqs - freq_mac)
colormap('jet');
colorbar()
xlabel('p')
ylabel('N')
title('Bursting Frequency')
xticks(1:length(Ps))
xticklabels(Ps)
yticks(1:length(Ns))
yticklabels(Ns)

figure()        
imagesc(amps - amp_mac)
colormap('jet');
colorbar()
xlabel('p')
ylabel('N')
title('Bursting Amplitude')
xticks(1:length(Ps))
xticklabels(Ps)
yticks(1:length(Ns))
yticklabels(Ns)

% END OF FILE