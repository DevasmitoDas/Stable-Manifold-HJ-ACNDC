%%%% 緩和率の感度 / hysteresis relaxation sweep (paper Remark 2, Fig. 3) %%%%
% For each lambda_h in cfg.relaxationRates the plant is changed, the problem is
% rebuilt (symbolic or analytic, buildProblem.m), the augmented LQR is RETUNED
% (f_design.m) and simulated from rest with B_d = 0.
% Remark 2: the slowest closed-loop pole is close to
%     -lambda_h/(1 + Ah*kh/bp)   (= -0.01167 for lambda_h = 0.02)
% Output : relaxCases{k}, relaxTable = [lambda_h, slowest pole, predicted pole, h(240)]
% Called from P5_solution_check.m

relaxCases = cell(1, numel(cfg.relaxationRates));
relaxTable = zeros(numel(cfg.relaxationRates), 4);
for k_rl = 1:numel(cfg.relaxationRates)
	par_rl = par; par_rl.lambda_h = cfg.relaxationRates(k_rl);
	fprintf('Relaxation sweep: lambda_h = %g (retuned augmented LQR)\n', par_rl.lambda_h);
	prob_rl = buildProblem(par_rl, cfg);
	des_rl = f_design(prob_rl, cfg);
	relaxCases{k_rl} = f_simulate_case('LQR', par_rl, des_rl, [], cfg, 0, 0, cfg.closedLoopTFinal);
	relaxCases{k_rl}.name = sprintf('LQR_lambda_%g', par_rl.lambda_h);
	relaxCases{k_rl}.hStar = des_rl.xStar(4);
	relaxTable(k_rl,:) = [par_rl.lambda_h, max(real(des_rl.poles)), ...
		-par_rl.lambda_h/(1 + par_rl.Ah*par_rl.kh/par_rl.bp), ...
		relaxCases{k_rl}.x(end,4) - des_rl.xStar(4)];
end
clear k_rl par_rl prob_rl des_rl
