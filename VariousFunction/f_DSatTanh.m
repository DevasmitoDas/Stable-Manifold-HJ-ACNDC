function dy = f_DSatTanh(u, umin, umax)
% ==================================================================
% dy = f_DSatTanh(u, umin, umax)
%
% F_DSATTANH  slope of the smooth saturation f_SatTanh
%     sigma'(u) = 1 - tanh((u - u_m)/u_h)^2,   0 < sigma'(u) <= 1
%   At the equilibrium of the paper, sigma'(u*) = 0.96.
% ==================================================================
um = (umin + umax)/2;  uh = (umax - umin)/2;
dy = 1 - tanh((u - um)/uh).^2;
end
