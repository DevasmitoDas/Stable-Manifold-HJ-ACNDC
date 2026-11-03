% ==========================================================================
% SP_value_fit: 価値関数の近似 / fit of the positive-definite value ansatz
%
%   Vhat(y) = (y'*P5*y)*exp(g(y)),  g = sum_{1<=|alpha|<=3} c_alpha y^alpha (55 coefficients)
% Collocation set (paper Sec. 4.1), all restricted to ||y|| <= 1.6:
%   - multi-scale Halton points (scales 0.1, 0.25, 0.5, 1 on [-0.8, 0.8]^5)
%   - 9 points on each coordinate axis
%   - 250 states of the augmented-LQR closed loop (cold start, B_d = 0)
% Objective eq. (fit): weighted HJB residual + manifold anchors + ridge,
% Gauss-Newton with frozen weights (polybase/f_value_fit.m).
% (Replaces SP_multiPoly of the reference toolbox; meshFlag = 0.)
% P5_solution_check.m で使用
% ==========================================================================

	% === LQR baseline (its closed-loop states enter the collocation set) ===
	if ~exist('caseLQR','var') || isempty(caseLQR)
		fprintf('\nSimulating the augmented LQR baseline (B_d = 0) ...\n');
		caseLQR = f_simulate_case('LQR', par, des, [], cfg, 0, 0, cfg.closedLoopTFinal);
		caseLQR.name = 'LQR_nominal';
	end
	if ~strcmp(caseLQR.status, 'completed')
		error('HJ:Baseline', 'Augmented LQR simulation failed.');
	end

	% === collocation set ===
	Ycoll = 1.6*f_halton(cfg.collocationPoints, [2 3 5 7 11], 37) - 0.8;
	pat_coll = [0.1 0.25 0.25 0.5 0.5 1 1 1 1 1];
	Ycoll = Ycoll.*pat_coll(1 + mod(0:cfg.collocationPoints-1, 10));
	ix_coll = unique(round(linspace(1, size(caseLQR.z,1), min(250, size(caseLQR.z,1)))));
	pathY = (caseLQR.z(ix_coll,:)' - des.zStar)./des.scale;
	for j_ax = 1:dim
		axisY = zeros(dim, 9); axisY(j_ax,:) = linspace(-0.8, 0.8, 9);
		Ycoll = [Ycoll axisY]; %#ok<AGROW>
	end
	Ycoll = [Ycoll pathY];
	Ycoll = Ycoll(:, sqrt(sum(Ycoll.^2,1)) <= cfg.hjDomainRadius);
	if isempty(Ycoll), error('HJ:FitDomain', 'Empty five-state collocation domain.'); end
	clear pat_coll ix_coll pathY axisY j_ax

	fprintf('\nFitting the five-variable scalar value function (%d collocation states) ...\n', size(Ycoll,2));
	valueModel = f_value_fit(des, manifold, Ycoll, cfg);

	if dbFlag == 1
		figure(71); clf
		semilogy(valueModel.fitHistory(:,1), valueModel.fitHistory(:,2), 'o-'); grid on
		xlabel('Gauss-Newton step'); ylabel('weighted HJB + anchor objective')
		title('Value fit (paper: 7363 -> 284)')
	end
