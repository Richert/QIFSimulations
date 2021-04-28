function [fr, spike_t, spike_n, mask, wait] = get_network_input(spike_t, wait, dt)
%GET_NETWORK_INPUT Summary of this function goes here
%   Detailed explanation goes here
spike_n = 1.*((spike_t>0)&(spike_t<=1));
mask = 1.*(round(wait)==0);
wait = max(wait-1,0);
spike_t = max(spike_t-1,0);
fr = spike_n ./ dt;
end