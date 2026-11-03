%%%% 周期外乱 / periodic disturbance sweep (paper Sec. 5.4, Table 3) %%%%
% w = 20*0.05*sin(2*pi*0.4*t) filtered by the lightly damped 0.7 Hz oscillator
% (zeta = 0.15); d = Cd*eta enters the angle equation through B_d/J.
% Both controllers face identical forcing for every B_d in cfg.disturbanceCouplings.
% Output : distCases{ib,ic} (ic = 1 LQR, 2 guarded HJ),
%          distTable = [B_d, J LQR, J HJ, dJ %, tail LQR, tail HJ, HJ use %, peak torque]
% Called from P5_solution_check.m

kinds_ds = {'LQR', 'HJ'};
distCases = cell(numel(cfg.disturbanceCouplings), 2);
distTable = zeros(numel(cfg.disturbanceCouplings), 8);
for ib_ds = 1:numel(cfg.disturbanceCouplings)
	bd_ds = cfg.disturbanceCouplings(ib_ds);
	for ic_ds = 1:2
		fprintf('Feedback %s: Bd=%g\n', kinds_ds{ic_ds}, bd_ds);
		distCases{ib_ds,ic_ds} = f_simulate_case(kinds_ds{ic_ds}, par, des, valueModel, cfg, ...
			bd_ds, cfg.disturbanceAmplitudeScale, cfg.closedLoopTFinal);
		distCases{ib_ds,ic_ds}.name = sprintf('%s_Bd_%g', kinds_ds{ic_ds}, bd_ds);
	end
	mL_ds = distCases{ib_ds,1}.metrics; mH_ds = distCases{ib_ds,2}.metrics;
	distTable(ib_ds,:) = [bd_ds, mL_ds.finiteHorizonCost, mH_ds.finiteHorizonCost, ...
		100*(mH_ds.finiteHorizonCost/mL_ds.finiteHorizonCost - 1), ...
		mL_ds.tailThetaRMSE, mH_ds.tailThetaRMSE, 100*mH_ds.hjUseFraction, mL_ds.peakDisturbanceTorque];
end
fprintf('\nTable 3: [Bd, J_LQR, J_HJ, dJ%%, tail LQR, tail HJ, HJ use %%, peak torque]\n'); disp(distTable)
clear kinds_ds ib_ds ic_ds bd_ds mL_ds mH_ds
