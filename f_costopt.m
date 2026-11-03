function Jcul = f_costopt(t, x, u, Q, R, Tf)
% ==================================================================
% Jcul = f_costopt(t, x, u, Q, R, Tf)
%
% F_COSTOPT  cost by the trapezoidal rule (reference toolbox routine)
%   J = int_0^Tf x'*Q*x + u'*R*u dt,  rows of x and u are samples.
% In this problem it is only a CROSS-CHECK: the costs of Tables 2-3 use the
% RK4-stage quadrature of the state cost and the exact held-rate energy
% (f_simulate_case.m), which is what the paper reports.
% ==================================================================
if nargin < 6, Tf = t(end); end
keep = t <= Tf;
t = t(keep); x = x(keep,:); u = u(keep,:);
L = sum((x*Q).*x, 2) + sum((u*R).*u, 2);
Jcul = trapz(t, L);
end
