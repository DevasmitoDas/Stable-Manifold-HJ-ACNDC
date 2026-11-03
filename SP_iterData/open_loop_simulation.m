%%%% open-loopのシミュレーション / open-loop response (paper Fig. 1) %%%%
% System : zdot = F_sigma(z) + G*nu + E*d  with the prescribed APPLIED command
%          sigma(u(t)) = 0.5 + 0.35*sin(pi*t/2) over 40 s and weak forcing B_d = 0.15
% The raw rate nu tracks sigma^{-1}(profile) at the samples (held-rate loop of
% f_simulate_case.m), so the applied command follows the profile exactly at t_k.
% Output : openCases{k} (one per cfg.openLoopProfiles), used for open_loop_loops.png
% Called from P5_solution_check.m

openCases = cell(1, numel(cfg.openLoopProfiles));
for k_ol = 1:numel(cfg.openLoopProfiles)
	fprintf('Open-loop response: %s (B_d = %g, %g s)\n', cfg.openLoopProfiles{k_ol}, par.Bd, cfg.openLoopTFinal);
	openCases{k_ol} = f_simulate_case(cfg.openLoopProfiles{k_ol}, par, des, [], cfg, par.Bd, 1, cfg.openLoopTFinal);
	openCases{k_ol}.name = ['open_loop_' cfg.openLoopProfiles{k_ol}];
end
clear k_ol
