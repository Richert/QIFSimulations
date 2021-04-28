close all; clear all;

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 5e-4;
dt2 = 5e-3;
T_sim = 410.0;
T = round(T_sim/dt);
T_init = 10.0;
T_rec = round((T_sim)/dt);
T_store = 0.0/dt;
v_mac_init = -2;

% QIF parameters
tau = 1.0;
Delta = 2.0*tau^2;
eta = -3.0*Delta;
J = 15.*sqrt(Delta);

% synpatic efficacy adaptation
tau_a = 10.0*tau;
alpha = 0.5/tau_a;

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
mu = 0.25*Delta; % -> for the low activity regime (eta = -6.0)
stim = 100+T_init;
dur = 200.0;
It = dt.*[1:T];
I = 0.*((It>0)&(It<T_sim)) + mu.*((It>stim)&(It<stim+dur));

for t = 1:T

    % calculate input
    IJ_mac = J.*r_mac_1.*(1 - min(1,e_mac_1));
    
    % macroscopic state variable evolution
    r_mac_2 = r_mac_1 + dt .* (Delta./(pi*tau) + 2 .* v_mac_1.*r_mac_1)/tau;
    v_mac_2 = v_mac_1 + dt .* (v_mac_1.^2 + eta + I(t) + IJ_mac*tau - (pi.*r_mac_1*tau).^2)/tau;   
    e_mac_2 = e_mac_1 + dt .* a_mac_1 / tau_a;
    a_mac_2 = a_mac_1 + dt .* (alpha*tau_a.*r_mac_1 - 2.*a_mac_1 - e_mac_1) / tau_a;
    
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

figure(); plot(r_mac_rec_av); title('r');
xticks([0, 100, 200, 300, 400])
xticklabels([0, 10, 20, 30, 40])
set(gca, 'PlotBoxAspectRatio',[4 1 1]);
figure(); plot(v_mac_rec_av); title('v');
xticks([0, 100, 200, 300, 400])
xticklabels([0, 10, 20, 30, 40])
set(gca, 'PlotBoxAspectRatio',[4 1 1]);
