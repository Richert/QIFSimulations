% Here we simulate the microscopic dynamics of the QIF model.

close all; clear all;
%dispstat('','init');
%dispstat('Start.','timestamp','keepthis')

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 1e-3;
dt2 = 1e-2;
N = 1e4;
N2 = 5e1;
idx = 1000;
T_sim = 40.0;
T = round(T_sim/dt);
T_init = 10.0;
T_rec = round(T_sim/dt);
T_store = 0.0;
v_mac_init = -2;
p = 1.0;

% QIF parameters
tau = 2.0;
Delta = 2.0*tau^2;
eta = -4.*Delta;
J = 15.*sqrt(Delta);
V_th = 1e2;

% eta adaptation 
eta_0 = eta+Delta.*tan((pi/2).*(2.*[1:N]-N-1)./(N+1));
%eta_0 = normrnd(eta,Delta,1,N);
tau_a = 10.0;
alpha = 0.05;

% initializations
%%%%%%%%%%%%%%%%%

v_1 = v_mac_init.*ones(size(eta_0)); % <-- optimal
e_1 = zeros(size(v_1));
a_1 = zeros(size(v_1));
wait = zeros(1,N);
spike_n = zeros(1,N);
spike_t = wait;

if p < 1
    spike_n = sparse(spike_n);
end
            
for m = 1:N; raster{m} = []; end

r_rec = zeros(1,T);
v_rec = zeros(size(r_rec));
e_rec = zeros(size(r_rec));
r_rec_av = zeros(1,T_rec.*dt2);
v_rec_av = zeros(size(r_rec_av));
e_rec_av = zeros(size(r_rec_av));
v_i_rec = zeros(size(r_rec));

% input definition
mu = 2.5*Delta; % -> for the low activity regime (eta = -6.0)
stim = T_init;
dur = 20.0;
It = dt.*[1:T];
I = 0.*((It>0)&(It<T_sim)) + mu.*((It>stim)&(It<stim+dur));

if p < 1.0
    C = rand(N, N);
    C_sorted = sort(C(:), 'ascend');
    C = sparse((C <= C_sorted(N^2*p)) * 1.0);
end

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
    
    % calculate state evolution of QIFs
    v_2 = v_1 + dt.*mask.*((v_1.^2 + eta_0 + I(t))./tau + IJ);
    e_2 = e_1 + dt .* mask .* a_1;
    a_2 = a_1 + dt .* mask .* (alpha.*s./tau_a - 2.*a_1./tau_a - e_1./tau_a.^2);
    
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
    v_i_rec(t) = v_2(idx);
    
    % coarse grained observations
    if t > T_store

        t_tmp = t-T_store;

        if (mod(t,1/dt2)==0)

            % microscopic state variables
            v_rec_av(t_tmp.*dt2) = mean(v_rec(t-1/dt2+1:t));
            r_rec_av(t_tmp.*dt2) = mean(r_rec(t-1/dt2+1:t));
            e_rec_av(t_tmp.*dt2) = mean(e_rec(t-1/dt2+1:t));

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

figure; plot(v_rec_av); title('v');
xticks([0, 100, 200, 300, 400])
xticklabels([0, 10, 20, 30, 40])
set(gca, 'PlotBoxAspectRatio',[4 1 1]);
figure; plot(r_rec_av); title('r');
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

% END OF FILE