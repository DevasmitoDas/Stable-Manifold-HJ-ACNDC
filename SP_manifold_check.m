% ==========================================================================
% SP_manifold_check: acceptance tests of one manifold trajectory (paper Table 1)
%
%   relH      relative |Hbar| along the trajectory and at the interval midpoints
%             |Hbar| / (l + |mu'*(a_sigma + G_y nu*)|)             < 5e-3
%   defect    midpoint forcing defect (linear interpolation of n)  < 5e-7
%   fpres     fixed-point relative residual (from SP_appNth)       < 1e-7 (5*tol)
%   tailNorm  terminal norm ||y(T_m)||                             < max(1e-7, 0.01*radius)
%   boundary  ||s(0) - xi|| + ||chi(T_m)||                         < 1e-10
%
% Uses : z, ini, radius, ker, des, cfg, dim      Sets : accepted, relH, defect,
%        tailNorm, boundary
% P4_iterative_calcu.m で使用
% ==========================================================================

src_k = i_nonlinear_term(z, des);
mid = f_manifold_midpoints(z, src_k, ker);
forcingMid = 0.5*(src_k(:,1:end-1) + src_k(:,2:end));
defect = max(sqrt(sum((i_nonlinear_term(mid, des) - forcingMid).^2, 1)));
yp_chk = Diagonalize*[z mid];
[H_chk, stage_chk, term_chk] = f_hamcalcu(yp_chk(1:dim,:), yp_chk(dim+1:end,:), des);
relH = abs(H_chk)./(1e-12 + stage_chk + abs(term_chk));
terminal = Diagonalize*z(:,end);
tailNorm = norm(terminal(1:dim));
boundary = norm(z(1:dim,1) - ini) + norm(z(dim+1:end,end));
accepted = all(isfinite(yp_chk(:))) && max(relH) < cfg.manifoldHamiltonianTolerance ...
	&& defect < cfg.manifoldDefectTolerance && boundary < 1e-10 ...
	&& tailNorm < max(1e-7, 0.01*radius);
clear src_k mid forcingMid yp_chk H_chk stage_chk term_chk terminal
