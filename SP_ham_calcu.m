% ==========================================================================
% SP_ham_calcu: Hamiltonian計算 / Hamiltonian on the manifold anchors
%
% HAMK : Hbar(y, mu) on every anchor (txpsumT), HAMKrel : relative residual
%        |Hbar|/(l + |mu'*(a_sigma + G_y nu*)|) (the valid data have HAMKrel << 1)
% P5_solution_check.m で使用
% ==========================================================================

	[HAMK, HAMKstage, HAMKterm] = f_hamcalcu(manifold.Y, manifold.Mu, des);
	HAMKrel = abs(HAMK)./(1e-12 + HAMKstage + abs(HAMKterm));
	fprintf('Anchors: %d, max relative |Hbar| on anchors = %.3e\n', size(manifold.Y,2), max(HAMKrel));

	if dbFlag == 1
		figure(61); clf
		loglog(sqrt(sum(manifold.Y.^2,1)), max(HAMKrel,1e-18), 'r.')
		grid on
		xlabel('||y|| (scaled)')
		ylabel('relative |Hbar|')
		title('Hamiltonian on the stable-manifold anchors')
	end
