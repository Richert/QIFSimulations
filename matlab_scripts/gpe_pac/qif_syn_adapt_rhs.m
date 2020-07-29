function y_delta = qif_syn_adapt_rhs(ys,pars)
%QIF_SYN_ADAPT_RHS Calculates rhs of DDE system, representing two delay-coupled
%QIF populations.
%   One is excitatory (STN) and one inhibitory (GPe). They
%   are coupled via simple exponential, current-based synapses plus some
%   discrete, constant delay.
%   Input
%   -----
%   - ys: state matrix with 10 rows (number of state variables) and 3 columns (number of time-lags).
%   - pars: Vector with 16 system parameters, including the two time delays
%   (15, 16).

% STN update
y_delta(1,1) = (pars(16)/(pi*pars(9)) + 2.0*ys(1,1)*ys(2,1))/pars(9);
y_delta(2,1) = (ys(2,1)^2+pars(1))/pars(9) + ys(5,1) - ys(6,1) - pars(9)*(pi*ys(1,1))^2;

% GPe update
y_delta(3,1) = (pars(16)/(pi*pars(10)) + 2*ys(3,1)*ys(4,1))/pars(10); 
y_delta(4,1) = (ys(4,1)^2+pars(2))/pars(10) + ys(7,1) - ys(8,1) - ys(9,1) - pars(10)*(pi*ys(3,1))^2;

% synaptic updates
y_delta(5,1) = (pars(3)*pars(4)*ys(1,2)-ys(5,1))/pars(13);
y_delta(6,1) = (pars(3)*pars(5)*ys(3,3)-ys(6,1))/pars(14);
y_delta(7,1) = (pars(3)*pars(6)*ys(1,3)-ys(7,1))/pars(13);
y_delta(8,1) = (pars(3)*pars(7)*ys(3,2)-ys(8,1))/pars(14);

% spike-frequency adaptation
y_delta(9,1) = ys(10,1);
y_delta(10,1) = (pars(8)*ys(3,1) - 2.0*ys(10,1) - ys(9,1)/pars(15))/pars(15);

end

