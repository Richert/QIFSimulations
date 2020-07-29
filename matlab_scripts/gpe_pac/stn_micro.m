% Here we simulate the microscopic dynamics of the QIF model.

close all; clear all;
%dispstat('','init');
%dispstat('Start.','timestamp','keepthis')

% parameter settings
%%%%%%%%%%%%%%%%%%%%

% general parameters
dt = 5e-3;
dt2 = 1e-1;
N = 1500;
N2 = 1e3;
T_sim = 1000;
T = round(T_sim/dt);
v_mac_init = -2;

% QIF parameters
J = 10.0;
tau = 6.0;
V_th = 1e2;
deltas = 0.05:0.05:1.51;
etas = -10.0:1:-2.9;

% start simulations
mean_firing_rates = zeros(length(deltas), length(etas));
std_firing_rates = zeros(size(mean_firing_rates));
n = 1;
for i=1:length(deltas)
    for j=1:length(etas)
        
        Delta = deltas(i);
        eta = etas(j);
        eta_0 = eta+Delta.*tan((pi/2).*(2.*[1:N]-N-1)./(N+1));

        % initializations
        %%%%%%%%%%%%%%%%%

        v_1 = v_mac_init.*ones(size(eta_0)); % <-- optimal

        wait = zeros(1,N);
        spike_n = zeros(1,N);
        spike_t = wait;

        r_rec = zeros(1,T);
        spikes_rec = zeros(N,T);
        r_rec_av = zeros(1,T.*dt2);
        spikes_rec_av = zeros(N, T.*dt2);

        % simulation
        %%%%%%%%%%%%

        for t = 1:T

            % calculate current spikes
            wait = max(wait-1,0);
            mask = 1.*(wait==0);
            spike_n = 1.*(spike_t==1);
            s = sum(spike_n)/N./dt;

            % calculate state evolution of QIFs
            v_2 = v_1 + dt.*mask.*((v_1.^2 + eta_0) ./ tau + J.*s);

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
            spikes_rec(:,t) = spike_n;

            % coarse grained observations
            if mod(t,1/dt2)==0

                %if mod(t,1/dt)==0
                %    disp(sprintf('Time: %d of %d done.',t.*dt,T.*dt));
                %end

                % state variables
                spikes_rec_av(:,t.*dt2) = mean(spikes_rec(:, t-1/dt2+1:t),2)*1e3;
                r_rec_av(t.*dt2) = mean(r_rec(t-1/dt2+1:t))*1e3;

            end
        end

        mean_firing_rates(i,j) = mean(r_rec_av(1, 500:end));
        std_firing_rates(i,j) = mean(std(spikes_rec_av,0,1));
        fprintf('Sweep: %d of %d done.',n,length(etas)*length(deltas))
        n = n+1;
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
imagesc(mean_firing_rates - 20.0)
colorbar()
title('Mean(r) - 20.0')
ylabel('Delta')
xlabel('Eta')
set(gca, 'YTick', 1:length(deltas), 'YTickLabel', deltas);
set(gca, 'XTick', 1:length(etas), 'XTickLabel', etas);

figure()
imagesc(std_firing_rates - 5.0)
colorbar()
title('Std(r) - 5.0')
ylabel('Delta')
xlabel('Eta')
set(gca, 'YTick', 1:length(deltas), 'YTickLabel', deltas);
set(gca, 'XTick', 1:length(etas), 'XTickLabel', etas);

% END OF FILE