function ua = f_command_profile(t, kind, par, Tfinal)
% ==================================================================
% ua = f_command_profile(t, kind, par, Tfinal)
%
% Prescribed APPLIED command sigma(u(t)) for open-loop runs.
%   'sinusoidal'  : u0 + uamp*sin(omega*t)  = 0.5 + 0.35*sin(pi*t/2)  (paper Fig. 1)
%   'smooth_step' : u0 + uamp/2*(1 + tanh((t - stepTime)/stepTau))
%   'chirp'       : u0 + uamp*sin(2*pi*(f0*t + (f1-f0)/(2*Tfinal)*t^2))
% ==================================================================
switch kind
    case 'sinusoidal'
        ua = par.u0 + par.uamp*sin(par.omega*t);
    case 'smooth_step'
        ua = par.u0 + 0.5*par.uamp*(1 + tanh((t - par.stepTime)/par.stepTau));
    case 'chirp'
        slope = (par.chirpF1 - par.chirpF0)/Tfinal;
        ua = par.u0 + par.uamp*sin(2*pi*(par.chirpF0*t + 0.5*slope*t.^2));
    otherwise
        error('HJ:Profile','Unknown open-loop profile: %s', kind);
end
end
