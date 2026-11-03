function dstate = i_state_eq_for_control(t, state, nu, par, des, amp)
% ==================================================================
% dstate = i_state_eq_for_control(t, state, nu, par, des, amp)
%
% Right-hand side of the sampled simulation (held raw rate nu):
%   state = [z; eta],  z = [theta; v; p; h; u_raw],  eta = exosystem state
%   zdot  = F_sigma(z) + G*nu + E*d,   d = Cd*eta,  E = (Bd/J)*e2   (eq. (aug))
%   etadot = Aeta*eta + Beta*w,  w = amp*wamp*sin(2*pi*wfreq*t)
% The drift is the SAME function used in the HJ design (des.prob.f), so
% sigma(u) is evaluated at every Runge-Kutta stage.
% ==================================================================
eta = state(6:7); d = par.Cd*eta;
y = (state(1:5) - des.zStar)./des.scale;
fz = des.scale.*(des.prob.f(y) + des.B*nu);
fz(2) = fz(2) + par.Bd*d/par.J;
w = amp*par.wamp*sin(2*pi*par.wfreq*t);
dstate = [fz; par.Aeta*eta + par.Beta*w];
end
