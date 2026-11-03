% ==================================================================
% Decide_xi.m  (iniFlag = 1)
%
% Seeds xi for the stable-manifold iteration, in the DIAGONALISED coordinates
% s (the stable block of (s,chi) = T^{-1}(y,mu)), as in the reference toolbox.
% Paper Sec. 3.2: the seeds are +-e_i (i = 1..5) and eight Halton directions,
% all of radius radi = 2e-4  ->  18 seeds.
%
% Defines
%   xi_dir : unit directions, one per column (dim-by-18)
%   xi_sum : seeds radi*xi_dir, one per ROW (reference convention)
% ==================================================================

xi_dir = [eye(dim), -eye(dim)];
xi_extra = 2*f_halton(cfg.manifoldRandomDirections, [2 3 5 7 11], 23) - 1;
xi_extra = xi_extra./max(sqrt(sum(xi_extra.^2,1)), eps);
xi_dir = [xi_dir, xi_extra];
xi_sum = (radi*xi_dir)';
clear xi_extra
