% ==================================================================
% P2: 近似回数 k と収束半径の確認 / convergence of the successive approximation
%
% For every real eigen-direction v of F (columns of Trs) the seed xi = radi*v
% is iterated with SP_appNth (k = 0, 1, 2, ...) until the fixed-point residual
% is below cfg.manifoldIterationRelTol. The relative Hamiltonian at t = 0,
%     |Hbar(y(0),mu(0))| / (l + |mu'*(a_sigma + G_y nu*)|),
% and the fixed-point residual are recorded for every k. Inside the radius of
% convergence both decrease with k; if they grow, radi is too large.
% (The paper uses radi = 2e-4; all 18 seeds then converge in 5-8 iterations.)
%
% The mesh and the exact exponential weights computed here are reused in P4.
% ==================================================================

% === kernels on the non-uniform mesh (cached for P4) ===
kerCache = cell(1, cfg.manifoldMeshRefinements + 1);
kerCache{1} = f_manifold_kernels(des, cfg, 0);
ker = kerCache{1};

init_sum_P2 = Trs./sqrt(sum(Trs.^2,1));      % unit eigen-directions (columns)
HamkN_P2 = zeros(1, size(init_sum_P2,2));     % relative |Hbar(t=0)| at the last k
iter_P2 = zeros(1, size(init_sum_P2,2));
conv_P2 = false(1, size(init_sum_P2,2));
Hhist_P2 = cell(1, size(init_sum_P2,2));
Rhist_P2 = cell(1, size(init_sum_P2,2));

for i_P2 = 1:size(init_sum_P2,2)
	ini = radi*init_sum_P2(:,i_P2);
	normSeed = norm(ini);
	% === k = 0 : linear solution s_0(t) = e^{Ft} xi, chi_0 = 0 ===
	z = zeros(2*dim, numel(ker.t));
	for it_t = 1:numel(ker.t)
		z(1:dim,it_t) = ker.linear(:,:,it_t)*ini;
	end
	yp0 = Diagonalize*z(:,1);
	[H0, st0, tm0] = f_hamcalcu(yp0(1:dim), yp0(dim+1:end), des);
	Hk = abs(H0)/(1e-12 + st0 + abs(tm0)); Rk = nan;
	knum_i = 0; iterConverged = false; iterFailed = false; fpres = inf;
	while knum_i < knum && ~iterConverged && ~iterFailed
		knum_i = knum_i + 1;
		SP_appNth   % N次近似計算 / one successive approximation step
		yp0 = Diagonalize*z(:,1);
		[H0, st0, tm0] = f_hamcalcu(yp0(1:dim), yp0(dim+1:end), des);
		Hk(end+1) = abs(H0)/(1e-12 + st0 + abs(tm0)); %#ok<SAGROW>
		Rk(end+1) = fpres; %#ok<SAGROW>
	end
	HamkN_P2(i_P2) = Hk(end); iter_P2(i_P2) = knum_i; conv_P2(i_P2) = iterConverged;
	Hhist_P2{i_P2} = Hk; Rhist_P2{i_P2} = Rk;
end
clear z yp0 H0 st0 tm0 Hk Rk it_t

if dbFlag == 1
	figure(21); clf
	color_P2 = lines(size(init_sum_P2,2));
	for i_P2 = 1:size(init_sum_P2,2)
		subplot(2,1,1); semilogy(0:numel(Hhist_P2{i_P2})-1, Hhist_P2{i_P2}, '.-', 'Color', color_P2(i_P2,:)); hold on
		subplot(2,1,2); semilogy(0:numel(Rhist_P2{i_P2})-1, Rhist_P2{i_P2}, '.-', 'Color', color_P2(i_P2,:)); hold on
	end
	subplot(2,1,1); grid on; xlabel('k'); ylabel('relative |Hbar| @ t=0'); title('P2: eigen-directions of F')
	subplot(2,1,2); grid on; xlabel('k'); ylabel('fixed-point residual')
end

% ------------- report --------------- %
Report_P2 = {'min & max of relative |Hbar(t=0)|', min(HamkN_P2), max(HamkN_P2); ...
	'iterations to converge (min, max)', min(iter_P2), max(iter_P2); ...
	'all eigen-directions converged', all(conv_P2), {}; ...
	'radius of xi', radi, {}; 'max num of iterations', knum, {}};
fprintf('Reports on P1 & P2:\n')
disp(Report_P1), disp(Report_P2)
