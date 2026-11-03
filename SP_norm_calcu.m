% ==========================================================================
% SP_norm_calcu: ノルムとコスト / norms, Table 2 metrics and cost comparison
%
% norm_x_* : scaled augmented error ||y(t)||,  norm_u_* : |nu(t)|
% cost_0th / cost_Nth : J_240 of the augmented LQR / guarded HJ (state + rate),
%   state part by RK4-stage quadrature, rate part sum_k Ts*Rnu*nu_k^2 (paper Sec. 5)
% cost_*_trap : trapezoidal cross-check with the reference routine f_costopt
% P5_solution_check.m で使用
% ==========================================================================

	norm_x_0th = sqrt(sum(((caseLQR.z - des.zStar')./des.scale').^2, 2));
	norm_x_Nth = sqrt(sum(((caseHJ.z  - des.zStar')./des.scale').^2, 2));
	norm_u_0th = abs(caseLQR.nu);
	norm_u_Nth = abs(caseHJ.nu);

	if dbFlag == 1
		figure(51); clf
		subplot(2,1,1); plot(caseLQR.t, norm_x_0th, 'b--', caseHJ.t, norm_x_Nth, 'r'); grid on
		ylabel('||y||'); legend('augmented LQR', 'guarded HJ')
		subplot(2,1,2); plot(caseLQR.t(1:end-1), norm_u_0th, 'b--', caseHJ.t(1:end-1), norm_u_Nth, 'r'); grid on
		ylabel('|nu|'); xlabel('Time (s)')
	end

	% === コスト計算 / cost ===
	cost_0th = caseLQR.metrics.finiteHorizonCost;
	cost_Nth = caseHJ.metrics.finiteHorizonCost;
	cost_comparison = (cost_Nth/cost_0th)*100;
	cost_0th_trap = f_costopt(caseLQR.t, caseLQR.z - des.zStar', [caseLQR.nu; caseLQR.nu(end)], des.Qphysical, des.R);
	cost_Nth_trap = f_costopt(caseHJ.t,  caseHJ.z  - des.zStar', [caseHJ.nu;  caseHJ.nu(end)],  des.Qphysical, des.R);

	fprintf('\n')
	fprintf('%s%0.5g\n', 'cost value about LINEAR function (augmented LQR) : ', cost_0th);
	fprintf('%s%0.5g\n', 'cost value about NONLINEAR function (guarded HJ)  : ', cost_Nth);
	fprintf('%s%0.5g%s\n', 'Cost Comparison : ( Nonlinear / linear ) is  ', cost_comparison, ' %');
	fprintf('(trapezoidal cross-check: %0.5g / %0.5g)\n', cost_0th_trap, cost_Nth_trap);

	% === Table 2 of the paper ===
	mL = caseLQR.metrics; mH = caseHJ.metrics;
	Report_P5_nominal = { ...
		'J_240 (state + rate)',              [mL.finiteHorizonCost mH.finiteHorizonCost], 'paper: 11.340, 8.979'; ...
		'  state part',                      [mL.stateCost mH.stateCost],                 'paper: 11.071, 7.573'; ...
		'  rate part',                       [mL.rateCost mH.rateCost],                   'paper: 0.269, 1.406'; ...
		'angle RMSE, full horizon (rad)',    [mL.thetaRMSE mH.thetaRMSE],                 'paper: 0.01174, 0.01236'; ...
		'angle RMSE, last 30 s (rad)',       [mL.tailThetaRMSE mH.tailThetaRMSE],         'paper: 3.47e-4, 2.27e-4'; ...
		'settling time +-0.007 rad (s)',     [mL.settlingTime mH.settlingTime],           'paper: 1.0, 20.5'; ...
		'||e_z(240)||',                      [mL.finalStateErrorNorm mH.finalStateErrorNorm], 'paper: 0.0222, 0.0120'; ...
		'HJ use (intervals of 24000)',       [mL.hjUsedIntervals mH.hjUsedIntervals],     'paper: 0, 18382 (76.6%)'; ...
		'HJ fallbacks: domain, residual',    [mH.fallbackDomain mH.fallbackResidual],     'paper: 10, 5608'; ...
		'HJ fallbacks: value, decrease guard', [mH.fallbackValue mH.fallbackGuard],       'paper: 0, 0'; ...
		'raw u min, applied sigma(u) min (HJ)', [mH.rawCommandMin mH.controlMin],         'paper: -0.146, 0.070'; ...
		'p and h errors at 240 s (LQR)',     [mL.finalPressureError mL.finalHysteresisError], 'paper: -0.0108, 0.0159'; ...
		'p and h errors at 240 s (HJ)',      [mH.finalPressureError mH.finalHysteresisError], 'paper: -0.0058, 0.0086'; ...
		'angle part 80*int e_theta^2',       [mL.angleStateCost mH.angleStateCost],       'paper: 2.65, 2.93'; ...
		'remaining state cost',              [mL.stateCost-mL.angleStateCost mH.stateCost-mH.angleStateCost], 'paper: 8.43, 4.64'};
	fprintf('\nTable 2 (columns: augmented LQR, guarded HJ):\n'); disp(Report_P5_nominal)
	Report.P5_nominal = Report_P5_nominal;
	clear mL mH
