% Here we simulate the microscopic dynamics of the QIF model.

close all; clear all;
%dispstat('','init');
%dispstat('Start.','timestamp','keepthis')

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 1e-3;
dt2 = 1e-1;
N = 1500;
N2 = 1e3;
T_sim = 100;
T = round(T_sim/dt);
v_mac_init = -2;

% QIF parameters
tau = 14;
V_th = 1e1;
deltas = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6];
etas = -2:2:8;

% start simulations
mean_firing_rates = zeros(length(deltas), length(etas));
std_firing_rates = zeros(size(mean_firing_rates));
for i=1:length(deltas)
    for j=1:length(etas)
        
        Delta = deltas(i);
        eta = etas(j)*Delta;
        eta_0 = eta+Delta.*tan((pi/2).*(2.*[1:N]-N-1)./(N+1));
        J = -5.0*sqrt(Delta);

        % initializations
        %%%%%%%%%%%%%%%%%

        v_1 = v_mac_init.*ones(size(eta_0)); % <-- optimal

        wait = zeros(1,N);
        spike_n = zeros(1,N);
        spike_t = wait;

        r_rec = zeros(1,T);
        v_rec = zeros(size(r_rec));
        r_rec_av = zeros(1,T.*dt2);
        v_rec_av = zeros(size(r_rec_av));

        I = dt.*[1:T];
        I = 0.*((I>0)&(I<T_sim));

        % simulation
        %%%%%%%%%%%%

        for t = 1:T

            % calculate current spikes
            wait = max(wait-1,0);
            mask = 1.*(wait==0);
            spike_n = 1.*(spike_t==1);
            s = sum(spike_n)/N./dt;

            % calculate state evolution of QIFs
            v_2 = v_1 + dt.*mask.*((v_1.^2 + eta_0 + I(t)) ./ tau + J.*s);

            % spiking mechanism
            spike = (v_2>V_th);
            wait(spike) = max(1, round((2.*tau./v_2(spike))./dt));
            spike_t(spike) = max(1, round((1.*tau./v_2(spike))./dt));
            spike_t = max(spike_t-1,0);
            v_2(spike) = -v_2(spike);

            % state variable updates
            v_1 = v_2;

            % variable recordings
            r_rec(t) = mean(spike_n)./dt;
            v_rec(t) = mean(v_2(wait==0));

            % coarse grained observations
            if mod(t,1/dt2)==0

                if mod(t,1/dt)==0
                    %dispstat(sprintf('Time: %d of %d done.',t.*dt,T.*dt),'timestamp');
                end

                % state variables
                v_rec_av(t.*dt2) = mean(v_rec(t-1/dt2+1:t));
                r_rec_av(t.*dt2) = mean(r_rec(t-1/dt2+1:t))*1e3;

            end
        end

        mean_firing_rates(i,j) = mean(r_rec_av(1, 500:end));
        std_firing_rates(i,j) = std(r_rec_av(1, 500:end));
        
    end
end

% plotting 
%%%%%%%%%%

figure()
imagesc(mean_firing_rates)
colorbar()
title('Mean(r)')
ylabel('Delta')
xlabel('Eta')
set(gca, 'YTick', 1:length(deltas), 'YTickLabel', deltas);
set(gca, 'XTick', 1:length(etas), 'XTickLabel', etas);

figure()
imagesc(std_firing_rates)
colorbar()
title('Std(r)')
ylabel('Delta')
xlabel('Eta')
set(gca, 'YTick', 1:length(deltas), 'YTickLabel', deltas);
set(gca, 'XTick', 1:length(etas), 'XTickLabel', etas);


figure()
imagesc((66.0 - mean_firing_rates).^2)
colorbar()
title('Mean(r) SE')
ylabel('Delta')
xlabel('Eta')
set(gca, 'YTick', 1:length(deltas), 'YTickLabel', deltas);
set(gca, 'XTick', 1:length(etas), 'XTickLabel', etas);

figure()
imagesc((19.0 - std_firing_rates).^2)
colorbar()
title('Std(r) SE')
ylabel('Delta')
xlabel('Eta')
set(gca, 'YTick', 1:length(deltas), 'YTickLabel', deltas);
set(gca, 'XTick', 1:length(etas), 'XTickLabel', etas);

% END OF FILE