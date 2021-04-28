function [v, spike_t, wait] = spiking_mechanism(v, v_th, IJ, etas, dt, wait, spike_t, tau)
%SPIKING_MECHANISM Summary of this function goes here
%   Detailed explanation goes here
spike = (v>v_th);
net_in = etas + IJ;
if length(IJ) > 1
    net_in = net_in(spike);
end
wait_tmp = (2*tau./v(spike))./dt - (6*net_in./v(spike).^3)./dt;
spike_t_tmp = (tau./v(spike))./dt;
wait(spike) = max(wait_tmp, spike_t_tmp + 1);
spike_t(spike) = spike_t_tmp;
v(spike) = -v(spike);
end

