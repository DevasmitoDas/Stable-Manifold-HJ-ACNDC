% ==========================================================================
% SP_response_calcu: 最適入力に対する応答の計算 / nominal closed loops (B_d = 0)
%
% Both controllers start from rest, z(0) = 0 (sigma(0) = 0.1192, ||y(0)|| = 1.733),
% Ts = 0.01 s for 240 s, held rate, two RK4 substeps per sample (paper Sec. 5).
%   caseLQR : augmented LQR  nu_L = sat_20(-K5*y)          (k = 0, "linear")
%   caseHJ  : guarded HJ selector of f_feedback_rate.m     ("nonlinear")
% P5_solution_check.m で使用
% ==========================================================================

	if ~exist('caseLQR','var') || isempty(caseLQR)
		caseLQR = f_simulate_case('LQR', par, des, [], cfg, 0, 0, cfg.closedLoopTFinal);
	end
	caseLQR.name = 'LQR_nominal';
	fprintf('\nSimulating the guarded HJ controller (B_d = 0) ...\n');
	caseHJ = f_simulate_case('HJ', par, des, valueModel, cfg, 0, 0, cfg.closedLoopTFinal);
	caseHJ.name = 'HJ_nominal';
	if ~strcmp(caseHJ.status, 'completed')
		warning('HJ:Response', 'Guarded HJ simulation diverged.');
	end

	if dbFlag == 1
		figure(41); clf
		subplot(2,1,1); plot(caseHJ.t, caseHJ.x); grid on
		title('Time response ( guarded HJ ): theta, v, p, h')
		subplot(2,1,2); stairs(caseHJ.t(1:end-1), caseHJ.nu); grid on; ylabel('nu'); xlabel('Time (s)')
		figure(42); clf
		subplot(2,1,1); plot(caseLQR.t, caseLQR.x); grid on
		title('Time response ( augmented LQR ): theta, v, p, h')
		subplot(2,1,2); stairs(caseLQR.t(1:end-1), caseLQR.nu); grid on; ylabel('nu'); xlabel('Time (s)')
	end
