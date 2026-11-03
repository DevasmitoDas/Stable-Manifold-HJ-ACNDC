% ==========================================================================
% SP_report: 数値のまとめ / numerical summary with the values of the paper
%
%   result/numerical_summary.txt   Reports P1-P5, Tables 1-3, relaxation sweep
%   result/closed_loop_metrics.csv nominal and disturbance cases
%   result/relaxation_metrics.csv, result/open_loop_metrics.csv
% P5_solution_check.m で使用
% ==========================================================================

	summaryFile = fullfile(cfg.resultFolder, 'numerical_summary.txt');
	fid_rp = fopen(summaryFile, 'w');
	if fid_rp < 0
		warning('HJ:Report', 'Cannot write %s; printing to the command window.', summaryFile);
		fid_rp = 1;
	end
	fprintf(fid_rp, 'SATURATION-AWARE STABLE-MANIFOLD CONTROL OF A HYSTERETIC SOFT PNEUMATIC ACTUATOR\n');
	fprintf(fid_rp, 'Sakamoto & Das, ACNDC 2027 -- reproduction run of %s\n', datestr(now));
	fprintf(fid_rp, 'Problem folder: %s\n', cfg.problemFolder);
	fprintf(fid_rp, 'Problem route : %s (%s)\n', prob.route, prob.status);
	fprintf(fid_rp, 'z = [theta v p h RAW u], nu = du/dt, sigma(u) = 0.5 + 0.5*tanh((u-0.5)/0.5)\n');
	fprintf(fid_rp, 'Qz = diag(%s), Rnu = %g, S = diag(%s)\n', num2str(cfg.Qdiag'), cfg.Rnu, num2str(cfg.scale'));
	fprintf(fid_rp, 'u* (raw) = %.10f, sigma(u*) = %g, sigma''(u*) = %g\n', des.uRawStar, des.uAppliedStar, ...
		des.prob.saturationDerivative(des.uRawStar));
	fprintf(fid_rp, 'Closed-loop poles (paper: -6.620+-11.619i, -12.493, -4.158, -0.01190):\n');
	for i_rp = 1:numel(des.poles)
		fprintf(fid_rp, '   %.10g %+.10gi\n', real(des.poles(i_rp)), imag(des.poles(i_rp)));
	end
	if isfield(checks, 'preHamiltonianSubstitutionError')
		fprintf(fid_rp, ['Symbolic pipeline: stationarity %.2g, substitution %.2g, chain rule %.2g ' ...
			'(paper: 1.1e-16, 7.1e-15, 2.2e-16)\n'], checks.symbolicStationarityError, ...
			checks.preHamiltonianSubstitutionError, checks.saturationCostateChainRuleError);
	end

	f_print_report(fid_rp, 'P1: local design', Report.P1);
	f_print_report(fid_rp, 'P2: convergence along the eigen-directions of F', Report.P2);
	f_print_report(fid_rp, 'P3: seeds', Report.P3);
	f_print_report(fid_rp, 'P4: stable-manifold trajectories (paper Table 1)', Report.P4);
	f_print_report(fid_rp, 'P5: value fit (paper Sec. 5.2)', Report.P5_fit);
	f_print_report(fid_rp, 'P5: nominal regulation, columns [augmented LQR, guarded HJ] (paper Table 2)', Report.P5_nominal);
	fprintf(fid_rp, '  Cost comparison (HJ / LQR) = %.2f %%  (paper: -20.8 %%)\n', cost_comparison);

	if exist('distTable','var') && ~isempty(distTable)
		paperT3 = [0.15 12.98 11.98 -7.7 0.01153 0.01153 71.5; ...
		           0.5  46.93 46.71 -0.5 0.03846 0.03835 86.0; ...
		           1   165.20 161.78 -2.1 0.07694 0.07608 92.9; ...
		           2   643.38 595.89 -7.4 0.15405 0.14756 85.7];
		fprintf(fid_rp, '\nPeriodic disturbance (paper Table 3)\n');
		fprintf(fid_rp, '  %-6s %-10s %-10s %-8s %-9s %-9s %-7s | paper: J LQR, J HJ, dJ, tail LQR, tail HJ, HJ use\n', ...
			'Bd', 'J LQR', 'J HJ', 'dJ %', 'tail LQR', 'tail HJ', 'HJ use');
		for i_rp = 1:size(distTable,1)
			fprintf(fid_rp, '  %-6g %-10.2f %-10.2f %-8.1f %-9.5f %-9.5f %-7.1f |', distTable(i_rp,1:7));
			j_rp = find(abs(paperT3(:,1) - distTable(i_rp,1)) < 1e-12, 1);
			if ~isempty(j_rp), fprintf(fid_rp, ' %g, %g, %g, %g, %g, %g', paperT3(j_rp,2:7)); end
			fprintf(fid_rp, '\n');
		end
		fprintf(fid_rp, '  peak disturbance torque at the largest Bd: %.3f (paper: 0.172 at Bd = 2)\n', distTable(end,8));
	end

	if exist('relaxTable','var') && ~isempty(relaxTable)
		paperH240 = [0.005 0.140; 0.02 0.0159; 0.05 2.2e-4; 0.10 1.8e-7];
		fprintf(fid_rp, '\nHysteresis relaxation sweep (paper Remark 2, Fig. 3)\n');
		fprintf(fid_rp, '  %-9s %-14s %-24s %-12s | paper h(240)\n', 'lambda_h', 'slowest pole', '-lambda_h/(1+Ah*kh/bp)', 'h(240)-h*');
		for i_rp = 1:size(relaxTable,1)
			fprintf(fid_rp, '  %-9g %-14.5f %-24.5f %-12.4g |', relaxTable(i_rp,:));
			j_rp = find(abs(paperH240(:,1) - relaxTable(i_rp,1)) < 1e-12, 1);
			if ~isempty(j_rp), fprintf(fid_rp, ' %g', paperH240(j_rp,2)); end
			fprintf(fid_rp, '\n');
		end
	end

	fprintf(fid_rp, '\nNotes\n');
	fprintf(fid_rp, '- Raw u is a state of the design; the applied command sigma(u) is bounded by construction.\n');
	fprintf(fid_rp, '- The rate bound |nu| <= 20 and the HJ/LQR switching are handled outside the HJB design.\n');
	fprintf(fid_rp, '- The guard (eq. (guard)) is pointwise at the samples for d = 0; it does not certify switching or disturbances.\n');
	fprintf(fid_rp, '- Exact costates exist only near the origin (anchors ||y|| <= 2.63e-4); elsewhere Vhat is shaped by the HJB residual.\n');
	fprintf(fid_rp, '- The guarded selector alternates between HJ and LQR around rho_H = 0.05; small floating-point\n');
	fprintf(fid_rp, '  differences (MATLAB release, BLAS, symbolic vs analytic route) can change individual\n');
	fprintf(fid_rp, '  switching decisions, so the HJ-use count and the HJ cost may differ in the last digits.\n');
	if fid_rp ~= 1, fclose(fid_rp); end
	fprintf('Numerical summary written to %s\n', summaryFile);

	% === CSV tables ===
	f_write_metrics_csv(fullfile(cfg.resultFolder, 'closed_loop_metrics.csv'), [{caseLQR, caseHJ}, reshape(distCases', 1, [])]);
	if ~isempty(relaxCases), f_write_metrics_csv(fullfile(cfg.resultFolder, 'relaxation_metrics.csv'), relaxCases); end
	if ~isempty(openCases),  f_write_metrics_csv(fullfile(cfg.resultFolder, 'open_loop_metrics.csv'), openCases); end
	clear fid_rp i_rp j_rp paperT3 paperH240
