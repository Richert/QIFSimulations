clear all;
close all;

% general parameters
dt = 1e-4;
dt2 = 1e-3;
n_neurons = 2;
n_plot = 1e2;
T_sim = 120.0;
T = round(T_sim/dt);
T_init = 20.0;
T_rec = round(T_init/dt);
T_buffer = round(1/dt2);
T_store = T_rec;

% QIF parameters
tau = 1.0;
eta = -1.0;
J = 5.0;
v_th = 1e2;

% synpatic parameters
tau_a = 1.0*tau;
tau_s = 2.0*tau;
alpha = 2.0;
Delta = 0.1;
mu = -2.0*Delta;
K = 2.0.*sqrt(Delta);

% initializations
%%%%%%%%%%%%%%%%%

% microscopic
v_1 = -2.0.*ones(1,n_neurons); % <-- optimal
I_1 = 0.0;
u_1 = -1.0;
c_1 = 0.0;
wait = zeros(1,n_neurons);
spike_n = zeros(1,n_neurons);
spike_t = wait;
net_in = zeros(1,n_neurons);
ext_in = zeros(size(net_in));
            
for m = 1:n_neurons; raster{m} = []; end

v_rec = zeros(n_neurons,T);
u_rec = zeros(1,T);
c_rec = zeros(size(u_rec));
I_rec = zeros(size(u_rec));
v_rec_av = zeros(n_neurons,T_rec.*dt2);
u_rec_av = zeros(1,T_rec.*dt2);
c_rec_av = zeros(size(u_rec_av));
I_rec_av = zeros(size(u_rec_av));

% input definition
mus = [5.0, 5.0];
stims = [20.0, 60.0] + T_init;
dur = 20.0;
It = dt*[1:T];
I = 0.*((It>0)&(It<T));
for s=1:length(stims)
    I = I + mus(1,s).*((It>stims(s))&(It<stims(s)+dur));
end

% simulation
%%%%%%%%%%%%

for t = 1:T
    
    % calculate inputs to system
    [fr, spike_t, spike_n, mask, wait] = get_network_input(spike_t, wait, dt);
    net_in(1,2) = I_1;
    ext_in(1,1) = I(t);
    
    % calculate state evolution of QIFs and vesicles
    v_2 = v_1 + dt.*mask.*((v_1.^2 + eta + ext_in)./tau + net_in);
    I_2 = I_1 + dt.*(J*c_1-I_1/tau_s);
    c_2 = c_1 + dt.*(Delta/(tau_a*pi) + 2*c_1*u_1) / tau_a;
    u_2 = u_1 + dt.*(u_1^2 + mu + alpha*fr(1) - tau_a*K*u_1 - (pi*c_1*tau_a)^2) / tau_a;
    
    % microscopic spiking mechanism
    [v_2, spike_t, wait] = spiking_mechanism(v_2, v_th, net_in, eta, dt, wait, spike_t, tau);
    
    % microscopic state variable updates
    v_1 = v_2;
    I_1 = I_2;
    c_1 = c_2;
    u_1 = u_2;
    
    % state variable recordings
    v_rec(:,t) = v_2;
    I_rec(:,t) = I_2;
    u_rec(t) = u_2;
    c_rec(t) = c_2;
    
    % coarse grained observations
    if t > T_store

        t_tmp = t-T_store;

        if (mod(t,1/dt2)==0)

            v_rec_av(:,t_tmp.*dt2) = mean(v_rec(:,t-1/dt2+1:t),2);
            I_rec_av(t_tmp.*dt2) = mean(I_rec(t-1/dt2+1:t));
            u_rec_av(t_tmp.*dt2) = mean(u_rec(t-1/dt2+1:t));
            c_rec_av(t_tmp.*dt2) = mean(c_rec(t-1/dt2+1:t));

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

% settings
plot_size = [10 10 500 200];
font_size = 12;
label_size = 16;
times = dt/dt2:dt/dt2:T_sim-T_init;
timeticks = 1:100:length(times);
timeticklabels = times(timeticks);

figure('DefaultAxesFontSize', font_size, 'DefaultAxesFontName', 'Roboto', ...
    'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
hold on 
plot(v_rec_av');
ylabel('v', 'FontSize', label_size)
xlabel('time', 'FontSize', label_size)
legend('qif1', 'qif2')
hold off

figure('DefaultAxesFontSize', font_size, 'DefaultAxesFontName', 'Roboto', ...
    'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
hold on
plot(u_rec_av);
ylabel('u', 'FontSize', label_size)
xlabel('time', 'FontSize', label_size)
hold off

figure('DefaultAxesFontSize', font_size, 'DefaultAxesFontName', 'Roboto', ...
    'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
hold on
plot(c_rec_av);
ylabel('c', 'FontSize', label_size)
xlabel('time', 'FontSize', label_size)
hold off

figure('DefaultAxesFontSize', font_size, 'DefaultAxesFontName', 'Roboto', ...
    'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
hold on
plot(I_rec_av);
ylabel('I', 'FontSize', label_size)
xlabel('time', 'FontSize', label_size)
hold off

figure('DefaultAxesFontSize', font_size, 'DefaultAxesFontName', 'Roboto', ...
    'DefaultAxesFontWeight', 'bold', 'Position', plot_size);
hold on;
for m = 1:n_neurons
    if ~isempty(raster{m})
        plot(raster{m},m.*ones(numel(raster{m},1)),'k.', 'MarkerSize', 10)
    end
end
xlabel('time', 'FontSize', label_size)
ylabel('neuron #', 'FontSize', label_size)
hold off