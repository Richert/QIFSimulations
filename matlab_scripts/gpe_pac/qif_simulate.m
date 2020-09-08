% Here we simulate the microscopic dynamics of the QIF model.

close all; clear all;
%dispstat('','init');
%dispstat('Start.','timestamp','keepthis')

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 1e-6;
dt2 = 1e-3;
T_sim = 2.0;
T = round(T_sim/dt);
v_mac_init = -2;

% QIF parameters
eta_e = -4.0;
eta_i = -8.0;
k = 1.0;
k_ee = 14.18;
k_ei = 11.36;
k_ie = 47.71;
k_ii = 7.61;
alpha = 2.99;

tau_e = 0.006;
tau_i = 0.014;
tau_s1 = 0.006;
tau_s2 = 0.002;
tau_ampa = 0.002;
tau_gabaa = 0.004;
delta = 2.0;
tau_a = 0.2;

% initializations
%%%%%%%%%%%%%%%%%

% macroscopic state variables
Ve_1 = 0.0;
Vi_1 = 0.0;
Re_1 = 0.0;
Ri_1 = 0.0;
Iee_1 = 0.0;
Iei_1 = 0.0;
Iie_1 = 0.0;
Iii_1 = 0.0;
Ia_1 = 0.0;
Xa_1 = 0.0;

% state variable recordings
re_rec = zeros(1,T.*dt2);
ri_rec = zeros(1,T.*dt2);

% state variable buffers
Re_buffer = zeros(1,int32(tau_s1/dt));
Ri_buffer = zeros(1,int32(tau_s1/dt));
idx_1 = int32(tau_s1/dt);
idx_2 = int32(tau_s2/dt);

% simulation
%%%%%%%%%%%%

for t = 1:T
    
    % stn
    Re_2 = Re_1 + dt * (delta/(pi*tau_e^2) + 2*Re_1*Ve_1/tau_e);
    Ve_2 = Ve_1 + dt * ((Ve_1^2 + eta_e)/tau_e + Iee_1 - Iei_1 - tau_e*(pi*Re_1)^2);
    
    % gpe
    Ri_2 = Ri_1 + dt * (delta/(pi*tau_i^2) + 2*Ri_1*Vi_1/tau_i);
    Vi_2 = Vi_1 + dt * ((Vi_1^2 + eta_i)/tau_i + Iie_1 - Iii_1 - Ia_1 - tau_i*(pi*Ri_1)^2);
    
    % synapses
    Iee_2 = Iee_1 + dt .* (k*k_ee*Re_buffer(idx_2) - Iee_1) / tau_ampa;
    Iei_2 = Iei_1 + dt .* (k*k_ei*Ri_buffer(idx_1) - Iei_1) / tau_gabaa;
    Iie_2 = Iie_1 + dt .* (k_ie*Re_buffer(idx_1) - Iie_1) / tau_ampa;
    Iii_2 = Iii_1 + dt .* (k*k_ii*Ri_buffer(idx_2) - Iii_1) / tau_gabaa;
    
    % sfa
    Ia_2 = Ia_1 + dt * Xa_1;
    Xa_2 = Xa_1 + dt * (alpha*Ri_1/tau_a - 2.0*Xa_1/tau_a - Ia_1/tau_a^2);
    
    % buffer update
    Re_buffer = circshift(Re_buffer,1);
    Ri_buffer = circshift(Ri_buffer,1);
    Re_buffer(1,1) = Re_2;
    Ri_buffer(1,1) = Ri_2;
    
    % macroscopic state variable updates
    Ve_1 = Ve_2;
    Re_1 = Re_2;
    Vi_1 = Vi_2;
    Ri_1 = Ri_2;
    Iee_1 = Iee_2;
    Iii_1 = Iii_2;
    Iei_1 = Iei_2;
    Iie_1 = Iie_2;
    Ia_1 = Ia_2;
    Xa_1 = Xa_2;
    
    % coarse grained observations
    if mod(t,1/dt2)==0
        
        % macroscopic state variables
        re_rec(t.*dt2) = Re_2;
        ri_rec(t.*dt2) = Ri_2;
        
    end
   
end


% plotting
%%%%%%%%%%

% firing rate
figure('DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Roboto', 'DefaultAxesFontWeight', 'bold');
hold on; 
plot(re_rec, 'LineWidth', 1.5, 'color', 'b');
plot(ri_rec, 'LineWidth', 1.5, 'color', 'r');
pbaspect([3 1 1])
xticks(linspace(1, T*dt2, 5));
xlabel('Time (s)')
ylabel('Firing rate (r)', 'FontSize', 16);
hold off;
