close all; clear all;

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 5e-4;
dt2 = 5e-3;
N2 = 5e1;
T_sim = 540.0;
T = round(T_sim/dt);
T_init = 20.0;
T_rec = round((T_sim - T_init)/dt);
T_store = T_init/dt;
v_mac_init = -2;

% QIF parameters
Delta = 2.0;
eta = -8.0;
J = 15.*sqrt(Delta);
V_th = 1e2;

% synpatic efficacy adaptation
tau_a = 10.0;
alpha = 0.05;

% sweep parameters
Ns = [1e3, 4e3];
Ps = [0.01, 0.1, 1.0];
reps = 1;
freqs = zeros(length(Ns), length(Ps));
amps = zeros(size(freqs));
fs = 1000;

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
mu = 2.5; % -> for the low activity regime (eta = -6.0)
stim = T_init;
dur = 500.0;
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

peaks_mac = findpeaks(r_mac_rec_av,'MinPeakProminence',0.5,'MinPeakDistance',100);
freq_mac = length(peaks_mac)/5;
amp_mac = max(r_mac_rec_av(1000:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% microscopic simulations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r=1:reps
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
            if p < 1
                spike_n = spike_n;
            end
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
                wait = max(wait-1,0);
                mask = 1.*(round(wait)==0);
                spike_n = 1.*((spike_t>0)&(spike_t<=1));
                spike_t = max(spike_t-1,0);
                if p < 1
                    s = spike_n*C./(dt*N);
                else
                    s = sum(spike_n)./(dt*N);
                end
                IJ = J.*s.*(1 - min(1, e_1))./p;

                % microscopic state variable evolution
                v_2 = v_1 + dt .* mask .* (v_1.^2 + eta_0 + I(t) + IJ);
                e_2 = e_1 + dt .* mask .* a_1;
                a_2 = a_1 + dt .* mask .* (alpha.*s./(tau_a*p) - 2.*a_1./tau_a - e_1./tau_a.^2);

                % microscopic spiking mechanism
                spike = (v_2>V_th);
                wait(spike) = (2./v_2(spike))./dt - (6*(eta_0(spike)+IJ(spike))./v_2(spike).^3)./dt;
                spike_t(spike) = (1./v_2(spike))./dt;
                v_2(spike) = -v_2(spike);

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

            peaks = findpeaks(r_rec_av,'MinPeakProminence',0.5,'MinPeakDistance',100);
            freqs(i,j) = length(peaks)/5;
            amps(i,j) = max(r_rec_av(1000:end));
            
            %figure()
            %hold on
            %plot(r_mac_rec_av(2000:4000))
            %plot(r_rec_av(2000:4000))
            %xlim([0,2000])
            %ylabel('r')
            %title(sprintf('N = %i, p = %.2f', N, p))
            %legend('macro', 'micro')
            %set(gca, 'PlotBoxAspectRatio',[4 1 1]);
            %hold off
            
            j = j+1;
            
        end
        i = i+1;
    end
end

set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000');
                      
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
saveas(gcf, 'burst_freq', 'svg')

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
saveas(gcf, 'burst_amp', 'svg')

% END OF FILE