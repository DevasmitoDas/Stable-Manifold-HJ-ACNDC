% ==================================================================
% P4: 与えたれたξに対して (y(t), mu(t)) を求める / stable-manifold trajectories
%
% For each seed xi (P3):
%   1. k = 0 : s_0(t) = e^{Ft} xi, chi_0 = 0 on the non-uniform mesh of [0, T_m]
%   2. successive approximations SP_appNth until the fixed point is reached
%   3. acceptance tests SP_manifold_check (paper Table 1)
%   4. if a test fails: refine the mesh once, then halve the radius (<= 4 times)
%   5. accepted trajectories are mapped back, (y,mu) = T*(s,chi), and sampled at
%      32 linear times on [0,2] and 32 logarithmic times on [2,T_m]:
%      these (y, mu) pairs are the manifold ANCHORS used by the value fit (P5).
% No negative-time extension is used (the anchors stay within ||y|| <= 2.63e-4).
%
% Outputs:
%   manifold   .Y, .Mu (anchors, columns), .t, .seedIndex and per-seed diagnostics
%   txpsumT    rows [t, y', mu'] of all anchors (reference "shooting lines")
%   txp_seedT  the same per seed (cell),   ini_sumT : (y(0), mu(0)) per seed
%   Report_P4  Table 1 of the paper
% ==================================================================

if ~exist('kerCache','var') || isempty(kerCache{1})
	kerCache = cell(1, cfg.manifoldMeshRefinements + 1);
	kerCache{1} = f_manifold_kernels(des, cfg, 0);
end

manifold = struct('Y',zeros(dim,0),'Mu',zeros(dim,0),'t',zeros(1,0),'seedIndex',zeros(1,0), ...
	'accepted',0,'attempted',ini_row,'seedRadius',[],'relativeHamiltonian',[], ...
	'terminalNorm',[],'defect',[],'iterations',[],'fixedPointRelativeResidual',[], ...
	'boundaryResidual',[],'horizon',des.horizon,'messages',{{}}, ...
	'method','Five-state Lyapunov-Perron stable forward/backward integral iteration');
txpsumT  = zeros(0, 1 + 2*dim);
txp_seedT = cell(1, ini_row);
ini_sumT = zeros(0, 2*dim);
NN = 0;

for xi_i = 1:ini_row
	accepted = false; message = 'No converged iteration.'; knum_i = 0; radius = radi;
	for retry = 1:cfg.manifoldRadiusRetries
		radius = radi/2^(retry - 1);
		ini = radius*xi_dir(:,xi_i);
		normSeed = norm(ini);
		for refinement = 0:cfg.manifoldMeshRefinements
			if isempty(kerCache{refinement + 1})
				kerCache{refinement + 1} = f_manifold_kernels(des, cfg, refinement);
			end
			ker = kerCache{refinement + 1};
			try
				% === k = 0 : 0次近似 (linear solution) ===
				z = zeros(2*dim, numel(ker.t));
				for it_t = 1:numel(ker.t)
					z(1:dim,it_t) = ker.linear(:,:,it_t)*ini;
				end
				% === k = 1, 2, ... : N次近似 ===
				knum_i = 0; iterConverged = false; iterFailed = false; fpres = inf;
				while knum_i < knum && ~iterConverged && ~iterFailed
					knum_i = knum_i + 1;
					NN = NN + 1;
					SP_appNth
				end
				if ~iterConverged
					message = sprintf('Iteration residual=%g.', fpres);
					break
				end
				SP_manifold_check
				if accepted
					% === anchors: 32 linear + 32 logarithmic sampling times ===
					desired = unique([linspace(0, min(2,des.horizon), 32), ...
						logspace(log10(min(2,des.horizon)), log10(des.horizon), 32)]);
					use = zeros(size(desired));
					for j_t = 1:numel(desired)
						[~, use(j_t)] = min(abs(ker.t - desired(j_t)));
					end
					use = unique(use);
					data = Diagonalize*z(:,use);                % (y, mu) = T*(s, chi)
					keep = sqrt(sum(data(1:dim,:).^2,1)) > 1e-8;
					data = data(:,keep); tdata = ker.t(use); tdata = tdata(keep);
					manifold.Y  = [manifold.Y  data(1:dim,:)];
					manifold.Mu = [manifold.Mu data(dim+1:end,:)];
					manifold.t  = [manifold.t tdata];
					manifold.seedIndex = [manifold.seedIndex xi_i*ones(1,size(data,2))];
					manifold.accepted = manifold.accepted + 1;
					manifold.seedRadius(end+1) = radius;
					manifold.relativeHamiltonian(end+1) = max(relH);
					manifold.terminalNorm(end+1) = tailNorm;
					manifold.defect(end+1) = defect;
					manifold.iterations(end+1) = knum_i;
					manifold.fixedPointRelativeResidual(end+1) = fpres;
					manifold.boundaryResidual(end+1) = boundary;
					txp_seedT{xi_i} = [tdata(:), data'];
					txpsumT = [txpsumT; txp_seedT{xi_i}]; %#ok<AGROW>
					y0mu0 = Diagonalize*z(:,1);
					ini_sumT = [ini_sumT; y0mu0']; %#ok<AGROW>
					message = 'accepted';
					break
				end
				message = sprintf('Hrel=%g, defect=%g, tail=%g.', max(relH), defect, tailNorm);
			catch ME_P4
				message = ME_P4.message;
				break
			end
		end
		if accepted, break; end
	end
	manifold.messages{end+1} = sprintf('Direction %d: %s', xi_i, message);
	fprintf('  Seed %d/%d: %s (radius %g, iterations %d)\n', xi_i, ini_row, message, radius, knum_i);
end
manifold.sampleRank = rank(manifold.Y, 1e-10);
clear z ker data keep tdata use desired y0mu0 j_t it_t

if manifold.accepted < cfg.manifoldMinAccepted || manifold.sampleRank < dim
	save(fullfile(cfg.resultFolder,'five_state_manifold_failure.mat'),'-v7');
	error('HJ:Manifold','Insufficient five-dimensional manifold data; failure checkpoint saved.');
end

% --------- report (paper Table 1) ---------------- %
Report_P4 = {'accepted seeds', manifold.accepted, sprintf('of %d (paper: 18 of 18)', ini_row); ...
	'relative |Hbar| along trajectory (worst)', max(manifold.relativeHamiltonian), 'threshold 5e-3, paper 1.28e-5'; ...
	'midpoint forcing defect (worst)', max(manifold.defect), 'threshold 5e-7, paper 2.13e-7'; ...
	'fixed-point relative residual (worst)', max(manifold.fixedPointRelativeResidual), 'threshold 1e-7, paper 1.34e-8'; ...
	'terminal norm ||y(T_m)|| (worst)', max(manifold.terminalNorm), 'threshold 2e-6, paper 7.24e-8'; ...
	'boundary residual (worst)', max(manifold.boundaryResidual), 'threshold 1e-10, paper 6.8e-21'; ...
	'iterations (min, max)', [min(manifold.iterations), max(manifold.iterations)], 'paper: 5-8'; ...
	'number of anchors', size(manifold.Y,2), 'paper: 1050'; ...
	'largest anchor radius ||y||', max(sqrt(sum(manifold.Y.^2,1))), 'paper: 2.63e-4'; ...
	'sample rank', manifold.sampleRank, 'must be 5'};
Report.P1 = Report_P1;
Report.P2 = Report_P2;
Report.P3 = Report_P3;
Report.P4 = Report_P4;
disp(Report.P3)
disp(Report.P4)
