% ===================================================================
% SP_model_evaluation: 近似性能の確認 / accuracy of the value fit (paper Sec. 5.2)
%
% 1. Relative gradient error at the manifold anchors
%        ||grad Vhat(y^m_j) - mu^m_j|| / max(1e-7, ||mu^m_j||)
% 2. HJB accuracy on held-out Halton states with ||y|| <= 1.6, measured by
%        rho_H = |R_H| / (1e-10 + l + |grad Vhat'*(a_sigma + G_y nu*)|)   (eq. (rho))
%    and the fraction that satisfies rho_H < 0.05 and |g| <= 3 (HJ eligible).
% 3. Structure at the origin: Vhat(0) = 0, grad Vhat(0) = 0 (Proposition 2).
% ===================================================================

	% === held-out states ===
	testY = 1.6*f_halton(4*cfg.validationPoints, [2 3 5 7 11], 10007) - 0.8;
	testY = testY(:, sqrt(sum(testY.^2,1)) <= cfg.hjDomainRadius);
	testY = testY(:, 1:min(cfg.validationPoints, size(testY,2)));
	if isempty(testY), error('HJ:ValidationDomain', 'Empty validation domain.'); end
	[~, ~, ~, Vtest, ~, gtest, rhoTest] = f_value_colloc(testY, valueModel, des, false);
	valueModel.validation = struct( ...
		'quantiles', [f_percentile(rhoTest,50) f_percentile(rhoTest,90) f_percentile(rhoTest,99)], ...
		'passFraction', mean(rhoTest < cfg.hjResidualMax & abs(gtest) <= cfg.maxLogCorrection), ...
		'passCount', sum(rhoTest < cfg.hjResidualMax & abs(gtest) <= cfg.maxLogCorrection), ...
		'residualPassFraction', mean(rhoTest < cfg.hjResidualMax), ...
		'minimumValue', min(Vtest), 'maximumLogCorrection', max(abs(gtest)), ...
		'points', size(testY,2), 'domainRadius', cfg.hjDomainRadius);

	% === anchor gradient errors ===
	[~, MuAnchor] = f_value_eval(manifold.Y, valueModel, des, false);
	gradErr = sqrt(sum((MuAnchor - manifold.Mu).^2,1))./max(1e-7, sqrt(sum(manifold.Mu.^2,1)));
	valueModel.anchorGradientError = [median(gradErr) f_percentile(gradErr,90) max(gradErr)];

	% === structure at the origin ===
	[V0, Mu0] = f_value_eval(zeros(dim,1), valueModel, des, false);
	valueModel.originCheck = [abs(V0) norm(Mu0)];

	% fitHistory rows: [step, objective before the step, accepted, line-search factor, ||delta||]
	nAccepted = sum(valueModel.fitHistory(:,3));
	idxLast = min(nAccepted + 1, size(valueModel.fitHistory,1));
	Report_P5_fit = {'anchor relative gradient error (median, p90, max)', valueModel.anchorGradientError, 'paper: 1.84e-5, 3.58e-5, 1.03e-4'; ...
		'objective (initial, after the accepted steps)', [valueModel.fitHistory(1,2), valueModel.fitHistory(idxLast,2)], 'paper: 7363 -> 284'; ...
		'accepted Gauss-Newton steps', nAccepted, 'paper: 9'; ...
		'collocation states', valueModel.trainingPoints, 'paper: 1294'; ...
		'held-out rho_H (median, p90, p99)', valueModel.validation.quantiles, 'paper: 0.099, 0.377, ~1'; ...
		'held-out HJ-eligible states', [valueModel.validation.passCount, valueModel.validation.points], 'paper: 153 of 500 (30.6%)'; ...
		'largest |g| on held-out states', valueModel.validation.maximumLogCorrection, 'paper: 3.19'; ...
		'Vhat(0), ||grad Vhat(0)||', valueModel.originCheck, 'must be 0'};
	disp(Report_P5_fit)
	Report.P5_fit = Report_P5_fit;
	clear MuAnchor V0 Mu0 nAccepted idxLast
