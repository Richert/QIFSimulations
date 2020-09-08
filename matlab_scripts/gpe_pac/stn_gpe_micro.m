% Here we simulate the microscopic dynamics of the QIF model.

close all; clear all;
%dispstat('','init');
%dispstat('Start.','timestamp','keepthis')

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 1e-2;
dt2 = 1e-1;
T_sim = 200;
T = round(T_sim/dt);
v_mac_init = -2;

% STN parameters
N_e = 15000;
tau_e = 6;
V_th_e = 1e2;
Delta_e = 0.4;
eta_e = -1.0;
eta_e = eta_e+Delta_e.*tan((pi/2).*(2.*[1:N_e]-N_e-1)./(N_e+1));
J_ee = 10.0;
J_ie = 100.0;

% GPe parameters
N_i = 45000;
tau_i = 14;
V_th_i = 1e2;
Delta_i = 0.4;
eta_i = 2.0;
eta_i = eta_i+Delta_i.*tan((pi/2).*(2.*[1:N_i]-N_i-1)./(N_i+1));
J_ii = 40.0;
J_ei = 80.0;

% connectivity matrices
C_ee = get_sparse_mat(N_e, N_e, 0.9, 1);
C_ie = get_sparse_mat(N_i, N_e, 0.9, 1);
C_ii = get_sparse_mat(N_i, N_i, 0.9, 1);
C_ei = get_sparse_mat(N_e, N_i, 0.9, 1);

% initializations
%%%%%%%%%%%%%%%%%

% STN
ve_1 = v_mac_init.*ones(size(eta_e));

wait_e = zeros(1,N_e);
spike_n_e = zeros(1,N_e);
spike_t_e = wait_e;

re_rec = zeros(1,T);
ve_rec = zeros(size(re_rec));
re_rec_av = zeros(1,T.*dt2);
ve_rec_av = zeros(size(re_rec_av));

% GPe
vi_1 = v_mac_init.*ones(size(eta_i));

wait_i = zeros(1,N_i);
spike_n_i = zeros(1,N_i);
spike_t_i = wait_i;

ri_rec = zeros(1,T);
vi_rec = zeros(size(ri_rec));
ri_rec_av = zeros(1,T.*dt2);
vi_rec_av = zeros(size(ri_rec_av));

% input
I = dt.*[1:T];
I = 0.*((I>0)&(I<T_sim));

% simulation
%%%%%%%%%%%%

for t = 1:T

    % calculate current spikes
    wait_e = max(wait_e-1,0);
    mask_e = 1.*(wait_e==0);
    spike_n_e = 1.*(spike_t_e==1);
    s_ee = (C_ee * (spike_n_e/N_e./dt)')';
    s_ie = (C_ie * (spike_n_e/N_e./dt)')';

    wait_i = max(wait_i-1,0);
    mask_i = 1.*(wait_i==0);
    spike_n_i = 1.*(spike_t_i==1);
    s_ei = (C_ei * (spike_n_i/N_i./dt)')';
    s_ii = (C_ii * (spike_n_i/N_i./dt)')';

    % calculate state evolution of QIFs
    ve_2 = ve_1 + dt.*mask_e.*((ve_1.^2 + eta_e + I(t)) ./ tau_e ...
        + J_ee .* s_ee - J_ei .* s_ei);
    vi_2 = vi_1 + dt.*mask_i.*((vi_1.^2 + eta_i) ./ tau_i ...
        + J_ie .* s_ie - J_ii .* s_ii);

    % spiking mechanism
    spike_e = (ve_2>V_th_e);
    wait_e(spike_e) = max(1, round((2.*tau_e./ve_2(spike_e))./dt));
    spike_t_e(spike_e) = max(1, round((1.*tau_e./ve_2(spike_e))./dt));
    spike_t_e = max(spike_t_e-1,0);
    ve_2(spike_e) = -ve_2(spike_e);

    spike_i = (vi_2>V_th_i);
    wait_i(spike_i) = max(1, round((2.*tau_i./vi_2(spike_i))./dt));
    spike_t_i(spike_i) = max(1, round((1.*tau_i./vi_2(spike_i))./dt));
    spike_t_i = max(spike_t_i-1,0);
    vi_2(spike_i) = -vi_2(spike_i);

    % state variable updates
    ve_1 = ve_2;
    vi_1 = vi_2;

    % variable recordings
    re_rec(t) = mean(spike_n_e)./dt;
    ve_rec(t) = mean(ve_2(wait_e==0));
    ri_rec(t) = mean(spike_n_i)./dt;
    vi_rec(t) = mean(vi_2(wait_i==0));

    % coarse grained observations
    if mod(t,1/dt2)==0
        ve_rec_av(t.*dt2) = mean(ve_rec(t-1/dt2+1:t));
        re_rec_av(t.*dt2) = mean(re_rec(t-1/dt2+1:t))*1e3;
        vi_rec_av(t.*dt2) = mean(vi_rec(t-1/dt2+1:t));
        ri_rec_av(t.*dt2) = mean(ri_rec(t-1/dt2+1:t))*1e3;
        if mod(t,1/dt)==0
            sprintf('Time: %d of %d done.',t.*dt,T.*dt)
        end
    end
end

% plotting 
%%%%%%%%%%

figure; plot(ve_rec_av); title('ve'); 
figure; plot(re_rec_av); title('re');
figure; plot(vi_rec_av); title('vi'); 
figure; plot(ri_rec_av); title('ri');

% END OF FILE