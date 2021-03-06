% Here we simulate the microscopic dynamics of the QIF model.

close all; clear all;
addpath('../')

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 1e-4;
dt2 = 1e-3;
N = 1e4;
N2 = 50;
T_sim = 450.0;
T = round(T_sim/dt);
T_init = 50.0;
T_rec = round((T_sim - T_init)/dt);
T_store = T_init/dt;
v_mac_init = -2;

% QIF parameters
tau = 1.0;
Delta = 2.0;
eta = 1.3;
J = 15.*sqrt(Delta);
eta_0 = eta+Delta.*tan((pi/2).*(2.*[1:N]-N-1)./(N+1));
v_th = 1000;

% adaptation 
tau_a = 10.0;
alpha = 1.0;

% initializations
%%%%%%%%%%%%%%%%%

% microscopic state variables
v_1 = v_mac_init.*ones(size(eta_0)); % <-- optimal
e_1 = zeros(size(v_1));
a_1 = zeros(size(v_1));
wait = zeros(1,N);
spike_n = zeros(1,N);
spike_t = wait;
            
for m = 1:N; raster{m} = []; end

% microscopic observation variables
r_rec = zeros(1,T);
v_rec = zeros(size(r_rec));
e_rec = zeros(size(r_rec));
r_rec_av = zeros(1,T_rec.*dt2);
v_rec_av = zeros(size(r_rec_av));
e_rec_av = zeros(size(r_rec_av));
v_i_rec = zeros(size(r_rec));

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
mus = [0.3, -1.9]; % -> for the bistable regime (eta = -4.6)
stims = [80., 240.] + T_init; % -> for the bistable regime (eta = -4.6)
% mus = [-1.5]; % -> for the high activity regime (eta = -4.0)
% stims = [100.0] + T_init; % -> for the high activity regime (eta = -4.0)
% mus = [0.5]; % -> for the low activity regime (eta = -6.0)
% stims = [100.0] + T_init; % -> for the low activity regime (eta = -6.0)
dur = 80.0;
It = dt.*[1:T];
I = 0.*((It>0)&(It<T_sim));
for s=1:length(stims)
    I = I + mus(1,s).*((It>stims(s))&(It<stims(s)+dur));
end

% simulation
%%%%%%%%%%%%

for t = 1:T
    
    % calculate current synaptic input
    [fr, spike_t, spike_n, mask, wait] = get_network_input(spike_t, wait, dt);
    s = mean(fr);
    
    % calculate state evolution of QIFs
    I_mic = eta_0 + I(t) + J.*s.*tau;
    v_2 = v_1 + dt.*mask.*((v_1.^2 + I_mic - e_1)./tau);
    e_2 = e_1 + dt .* a_1;
    a_2 = a_1 + dt .* (alpha.*fr./tau_a - 2.*a_1./tau_a - e_1./tau_a.^2);
    
    % macroscopic state variable evolution
    I_mac = eta + I(t) - e_mac_1 + J.*r_mac_1.*tau;
    r_mac_2 = r_mac_1 + dt .* (Delta./(pi*tau) + 2.*v_mac_1.*r_mac_1)/tau;
    v_mac_2 = v_mac_1 + dt .* (v_mac_1.^2 + I_mac - (pi.*r_mac_1*tau).^2)/tau;   
    e_mac_2 = e_mac_1 + dt .* a_mac_1 / tau_a;
    a_mac_2 = a_mac_1 + dt .* (alpha*tau_a.*r_mac_1 - 2.*a_mac_1 - e_mac_1) / tau_a;
    
    % microscopic spiking mechanism
    [v_2, spike_t, wait] = spiking_mechanism(v_2, v_th, dt, wait, spike_t, tau, I_mic);

    % microscopic state variable updates
    v_1 = v_2;
    a_1 = a_2;
    e_1 = e_2;
    
    % microscopic state variable recordings
    r_rec(t) = mean(spike_n)./dt;
    v_rec(t) = mean(v_2(wait==0));
    e_rec(t) = mean(e_2(wait==0));
%     v_i_rec(t) = v_2(idx);
    
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

            % microscopic state variables
            v_rec_av(t_tmp.*dt2) = mean(v_rec(t-1/dt2+1:t));
            r_rec_av(t_tmp.*dt2) = mean(r_rec(t-1/dt2+1:t));
            e_rec_av(t_tmp.*dt2) = mean(e_rec(t-1/dt2+1:t));
            
            % macroscopic state variables
            v_mac_rec_av(t_tmp.*dt2) = mean(v_mac_rec(t-1/dt2+1:t));
            r_mac_rec_av(t_tmp.*dt2) = mean(r_mac_rec(t-1/dt2+1:t));
            e_mac_rec_av(t_tmp.*dt2) = mean(e_mac_rec(t-1/dt2+1:t));
            
        end
        if sum(spike_n)>0
            rid = find(spike_n==1);
            for m = 1:numel(rid)
                raster{rid(m)} = cat(1,raster{rid(m)},t);
            end
        end
    end

end

% plotting
%%%%%%%%%%

figure; hold on; plot(v_rec_av); plot(v_mac_rec_av); hold off;
title('v');
xticks([0, 100, 200, 300, 400])
xticklabels([0, 10, 20, 30, 40])
set(gca, 'PlotBoxAspectRatio',[4 1 1]);
figure; hold on; plot(r_rec_av); plot(r_mac_rec_av); hold off;
title('r');
xticks([0, 100, 200, 300, 400])
xticklabels([0, 10, 20, 30, 40])
set(gca, 'PlotBoxAspectRatio',[4 1 1]);
figure; plot(v_i_rec); title('v_i');
xticks([0, 20000, 40000, 60000, 80000])
xticklabels([0, 10, 20, 30, 40])
set(gca, 'PlotBoxAspectRatio',[4 1 1]);

figure; hold on;
[~,randidx] = sort(randn(N,1));
for m = 1:N2
    id = randidx(m);
    if ~isempty(raster{id})
        plot(raster{id},m.*ones(numel(raster{id},1)),'k.')
    end
end
xticks([0, 20000, 40000, 60000, 80000])
xticklabels([0, 10, 20, 30, 40])
ylim([0, 51]);
set(gca, 'PlotBoxAspectRatio',[4 1 1]);

save('neco_fig5_data.mat', 'v_rec_av', 'r_rec_av', 'v_i_rec', 'raster', ...
    'v_mac_rec_av', 'r_mac_rec_av')

% END OF FILE