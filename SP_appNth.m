% ==========================================================================
% SP_appNth: N次近似計算スクリプト / one successive approximation k -> k+1
%
%   s_{k+1}(t)   = e^{Ft} xi + int_0^t e^{F(t-a)} n_s(s_k(a),chi_k(a)) da   (i_directforx)
%   chi_{k+1}(t) = -int_t^{T_m} e^{F'(a-t)} n_u(s_k(a),chi_k(a)) da         (i_directforp)
%   (paper eq. (iter), Sakamoto & van der Schaft 2008)
% The forcing is interpolated linearly on each mesh interval with exact
% exponential weights and the update is relaxed by cfg.manifoldRelaxation (0.8).
% Once the relative fixed-point residual is below tol, one more unrelaxed map
% is evaluated; the iteration has converged if that residual is below 5*tol.
%
% Uses : z (2n-by-nt, diagonalised coordinates [s; chi]), ini (seed xi), normSeed,
%        ker (f_manifold_kernels), des, cfg, dim
% Sets : z, fpres, iterConverged, iterFailed
% P2_hamiltonian_check.m, P4_iterative_calcu.m で使用
% ==========================================================================

src_k  = i_nonlinear_term(z, des);                         % [n_s; n_u] along iterate k
z_next = [i_directforx(ini, src_k(1:dim,:), ker); ...
          i_directforp(src_k(dim+1:2*dim,:), ker)];
if any(~isfinite(z_next(:))) || max(abs(z_next(:))) > 100*normSeed
	iterFailed = true;                                     % outside the radius of convergence
else
	fpres = max(abs(z_next(:) - z(:)))/max(normSeed, max(abs(z_next(:))));
	z = (1 - cfg.manifoldRelaxation)*z + cfg.manifoldRelaxation*z_next;
	if fpres < cfg.manifoldIterationRelTol
		src_k  = i_nonlinear_term(z, des);
		z_next = [i_directforx(ini, src_k(1:dim,:), ker); ...
		          i_directforp(src_k(dim+1:2*dim,:), ker)];
		fpres  = max(abs(z_next(:) - z(:)))/max(normSeed, max(abs(z_next(:))));
		iterConverged = fpres < 5*cfg.manifoldIterationRelTol;
	end
end
clear src_k z_next
