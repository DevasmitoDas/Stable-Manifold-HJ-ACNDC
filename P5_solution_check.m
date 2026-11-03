% ===================================================================
% P5: 解の確認 / value approximation, guarded feedback, experiments, figures
%
%  SP_ham_calcu         Hamiltonian on the manifold anchors
%  SP_shooting_plot     anchors projected on the (y_i, y_j) planes (dbFlag = 1)
%  SP_value_fit         collocation set + Gauss-Newton fit of Vhat (paper Sec. 4.1)
%  SP_model_evaluation  anchor gradient errors and held-out rho_H (paper Sec. 5.2)
%  SP_response_calcu    nominal augmented LQR vs guarded HJ (paper Sec. 5.3, Table 2)
%  SP_norm_calcu        norms, Table 2 metrics, cost comparison
%  SP_iterData/disturbance_sweep     Table 3
%  SP_iterData/relaxation_sweep      Fig. 3 (retuned LQR for each lambda_h)
%  SP_iterData/open_loop_simulation  Fig. 1
%  SP_graphplot_paper   the four figures of the paper (result/figures)
%  SP_report            numerical_summary.txt and CSV tables (with paper values)
% ===================================================================

	SP_ham_calcu
	if dbFlag == 1
		SP_shooting_plot
	end

	% === 価値関数の近似 / value function (meshFlag = 0) ===
	if meshFlag == 0
		SP_value_fit
	else
		error('HJ:P5', 'Only meshFlag = 0 (positive-definite polynomial-exponential ansatz) is implemented for this problem.');
	end
	SP_model_evaluation

	% === 最適制御入力の応答性 / closed-loop responses ===
	SP_response_calcu

	% === ノルム・コスト / norms and cost ===
	SP_norm_calcu

	% === 論文の実験 / experiments of the paper ===
	distCases = {}; relaxCases = {}; openCases = {};
	if cfg.runDisturbanceSweep
		disturbance_sweep
	end
	if cfg.runRelaxationSweep
		relaxation_sweep
	end
	if cfg.runOpenLoop
		open_loop_simulation
	end

	% === 図と数値のまとめ / figures and summary ===
	SP_graphplot_paper
	SP_report
