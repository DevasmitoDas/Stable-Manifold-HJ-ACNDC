function y = f_SatTanh(u, umin, umax)
% ==================================================================
% y = f_SatTanh(u, umin, umax)
%
% F_SATTANH  smooth (tanh) saturation of the actuator command (paper eq. (sat))
%     sigma(u) = u_m + u_h*tanh((u - u_m)/u_h),  u_m = (umin+umax)/2, u_h = (umax-umin)/2
%   sigma maps R onto (umin, umax) and 0 < sigma'(u) <= 1.
%   Unlike f_Sat / f_SatS of the reference toolbox, sigma models the actuator
%   itself and enters the drift of the augmented plant (it is not an input map).
%
% input  - u : raw command (any size),  umin, umax : applied-command range
% output - y : applied command sigma(u)
% ==================================================================
um = (umin + umax)/2;  uh = (umax - umin)/2;
y = um + uh*tanh((u - um)/uh);
end
