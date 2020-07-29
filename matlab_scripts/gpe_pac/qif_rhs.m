function y_delta = qif_rhs(ys,pars)
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
y_delta(1,1) = (pars(9)/(pi*pars(7)) + 2*ys(1,1)*ys(2,1))/pars(7);
y_delta(2,1) = (ys(2,1)^2+pars(1))/pars(7) + pars(3)*ys(1,2) - pars(4)*ys(3,3) - pars(7)*(pi*ys(1,1))^2;

% GPe update
y_delta(3,1) = (pars(10)/(pi*pars(8)) + 2*ys(3,1)*ys(4,1))/pars(8); 
y_delta(4,1) = (ys(4,1)^2+pars(2))/pars(8) + pars(5)*ys(1,3) - pars(6)*ys(3,2) - pars(8)*(pi*ys(3,1))^2;

end

