function [v, spike_t, wait] = spiking_mechanism(v, v_th, dt, wait, spike_t, tau, I)
%SPIKING_MECHANISM Summary of this function goes here
%   Detailed explanation goes here
spike = (v>v_th);
if any(spike)
    wait_tmp = (2*tau./v(spike))./dt - (6*I(spike)./v(spike).^3)./dt;
    spike_t_tmp = (tau./v(spike))./dt;
    wait(spike) = max(wait_tmp, 1);
    spike_t(spike) = max(spike_t_tmp, 1);
    v(spike) = -v(spike);
end
end