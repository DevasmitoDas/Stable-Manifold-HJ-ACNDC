% ===============================================================================
% 計算パラメータ設定 / Data settings  --  inputData/5D_Pneumatic_nu_sat_symb
%
% Saturation-aware stable-manifold optimal control of a hysteretic soft
% pneumatic actuator (N. Sakamoto, D. Das, ACNDC 2027, Documents/NS_DD_ACNDC_2027.pdf).
%
% Design state     z  = [theta; v; p; h; u]   (u = RAW command, integrator state)
% Optimised input  nu = du/dt                 (decision variable)
% Applied command  ua = sigma(u) = u_m + u_h*tanh((u-u_m)/u_h)   (part of the drift)
% Cost             J  = int (z-z*)'*Qz*(z-z*) + Rnu*nu^2 dt     (NO factor 1/2)
% Scaled coordinates y = S^{-1}(z-z*), scaled costate mu = S*lambda  (paper eq. (preH))
%
% This file is a SCRIPT (toolbox convention). It defines the option flags, the
% physical parameters (par), all design / simulation settings (cfg), the problem
% functions (prob) and the reference-toolbox names A, B, Q, R, dim, knum, radi, ...
%
% Problem functions (prob) are produced by buildProblem.m:
%   symbFlag = 1 : symbSetting.m -> symbCompile.m (Symbolic Math Toolbox), as in the paper
%   symbFlag = 0 : numericProblem.m (hand-derived, identical formulas, no toolbox)
% ===============================================================================

% ===============================================
% 各種フラグ設定 / Option flags (reference names)
% ===============================================
	satFlag  = 0;  % 0: the HJB problem contains no input-saturation function.
                   %    The actuator saturation sigma(u) is a STATE nonlinearity of the
                   %    drift (paper Sec. 2.2), so the problem stays input affine in nu.
	dbFlag   = 0;  % 1: show diagnostic figures (anchors, Hamiltonian, P2 convergence,
                   %    fit history, time responses). 0: only the four paper figures.
	meshFlag = 0;  % 0: polynomial representation. Here: the positive-definite ansatz
                   %    Vhat(y) = (y'*P5*y)*exp(g(y)) of polybase/ (paper eq. (ansatz)).
	kFlag    = 1;  % 1: maximal number of successive approximations knum from this file
	xiRFlag  = 1;  % 1: seed radius radi from this file
	xiTHFlag = 1;  % 1: no interactive angular division of the xi-sphere
	iniFlag  = 1;  % 1: the seeds xi are given in Decide_xi.m (+-e_i and Halton directions)
	symbFlag = 1;  % 1: symbolic derivation (Symbolic Math Toolbox), as in the paper
                   % 0: analytic formulas of numericProblem.m (automatic fallback)

% ===============================================
% 物理パラメータ / Actuator parameters (paper Sec. 2.1, normalised demonstration values)
% ===============================================
	par = struct('J',0.020,'k',0.600,'c',0.080,'bp',0.350, ...
	    'taup',0.200,'ku',1.0,'kh',0.250,'Ah',1.0,'beta',0.8, ...
	    'gamma',0.2,'n',1.0,'eps_abs_rate',1e-4,'eps_abs_h',1e-4, ...
	    'lambda_h',0.020,'hRelaxEq',0.0);

	% === Exosystem eta' = Aeta*eta + Beta*w, d = Cd*eta (used in simulation only) ===
	par.Bd    = 0.150;                 % torque coupling for the open-loop run (Fig. 1)
	par.wamp  = 0.050;                 % amplitude of w (multiplied by the experiment scale)
	par.wfreq = 0.40;                  % frequency of w [Hz]
	wd_eta = 2*pi*0.70; zeta_eta = 0.15;   % 0.7 Hz lightly damped oscillator
	par.Aeta = [0 1; -wd_eta^2 -2*zeta_eta*wd_eta];
	par.Beta = [0; 1];
	par.Cd   = [1 0];
	clear wd_eta zeta_eta

	% === Prescribed applied-command profiles for open-loop runs ===
	% sinusoidal: sigma(u) = u0 + uamp*sin(omega*t) = 0.5 + 0.35*sin(pi*t/2)  (Fig. 1)
	par.u0 = 0.5; par.uamp = 0.35; par.omega = 2*pi*0.25;
	par.stepTime = 5; par.stepTau = 0.4; par.chirpF0 = 0.05; par.chirpF1 = 0.80;

% ===============================================
% 制御問題設定 / Control problem (paper Sec. 3.1)
% ===============================================
	cfg = struct();
	cfg.thetaStar = 0.35;                  % target angle -> p* = 0.6, sigma(u*) = 0.6
	cfg.Qdiag     = [80; 5; 2; 1; 0.20];   % Qz = diag(80,5,2,1,0.2)   (no factor 1/2)
	cfg.Rnu       = 0.020;                 % rate weight R_nu
	cfg.scale     = [0.35; 1; 0.60; 0.40; 0.60];   % S = typical error sizes
	cfg.uMin      = 0;                     % sigma: R -> (uMin, uMax), u_m = u_h = 0.5
	cfg.uMax      = 1;
	cfg.uRateMax  = 20;                    % deployment bound |nu| <= nubar

% ===============================================
% 設計用パラメータ / Design parameters: stable manifold (paper Sec. 3.2)
% ===============================================
	cfg.manifoldRadius           = 2e-4;   % radius of the seeds xi
	cfg.manifoldRandomDirections = 8;      % Halton directions in addition to +-e_i
	cfg.manifoldMinAccepted      = 10;     % minimal number of accepted trajectories
	cfg.manifoldTailExponents    = 8;      % T_m = max(40, 8/|max Re eig(F)|)
	cfg.manifoldMinimumHorizon   = 40;
	cfg.manifoldFastStep         = 0.01;   % mesh steps on [0,4], [4,30], [30,T_m]
	cfg.manifoldMiddleStep       = 0.1;
	cfg.manifoldSlowStep         = 0.5;
	cfg.manifoldMaxIterations    = 80;     % maximal number of successive approximations
	cfg.manifoldIterationRelTol  = 1e-7;   % fixed-point relative residual
	cfg.manifoldRelaxation       = 0.8;    % relaxation of each update
	cfg.manifoldRadiusRetries    = 5;      % halve the radius at most 4 times
	cfg.manifoldMeshRefinements  = 1;      % halve the mesh at most once
	cfg.manifoldDefectTolerance  = 5e-7;   % midpoint forcing defect
	cfg.manifoldHamiltonianTolerance = 0.005;  % relative |Hbar| along the trajectory

% ===============================================
% 価値関数近似 / Value approximation (paper Sec. 4.1)
% ===============================================
	cfg.collocationPoints     = 1000;      % multi-scale Halton collocation points
	cfg.validationPoints      = 500;       % held-out Halton states
	cfg.valueCorrectionDegree = 3;         % g(y): 5 + 15 + 35 = 55 coefficients
	cfg.fitMaxIterations      = 25;        % Gauss-Newton steps
	cfg.fitRidge              = 0.005;     % lambda_c
	cfg.manifoldAnchorWeight  = 0.05;      % omega_a
	cfg.fitWeightFloor        = 0.03;      % w_i = 1/(0.03 + l_i)
	cfg.maxLogCorrection      = 3;         % |g| <= 3

% ===============================================
% ガード付き切替 / Guarded rate selection (paper Sec. 4.2)
% ===============================================
	cfg.hjDomainRadius     = 1.6;          % HJ eligible only if ||y|| <= 1.6
	cfg.hjResidualMax      = 0.05;         % rho_H <= 0.05
	cfg.hjDecreaseFraction = 0.02;         % decrease guard (eq. (guard))

% ===============================================
% 応答性検証用パラメータ / Simulation (paper Sec. 5)
% ===============================================
	cfg.Ts                = 0.01;          % held-rate sample time
	cfg.rkSubsteps        = 2;             % RK4 substeps per sample
	cfg.closedLoopTFinal  = 240;           % [s]
	cfg.openLoopTFinal    = 40;            % [s] (Fig. 1)
	cfg.initialState      = [];            % [] -> z(0) = 0 (cold start from rest)
	cfg.initialRawCommand = 0;             % raw u(0) = 0, sigma(0) = 0.1192
	cfg.initialEta        = [0; 0];
	cfg.tailWindow        = 30;            % tail RMSE over the last 30 s
	cfg.settlingBand      = 0.007;         % |theta - theta*| <= 0.007 rad
	cfg.divergenceScaledNorm = 30;

	% === Experiments ===
	cfg.runDisturbanceSweep       = true;  % Table 3
	cfg.disturbanceCouplings      = [0.15 0.5 1 2];
	cfg.disturbanceAmplitudeScale = 20;
	cfg.runRelaxationSweep        = true;  % Fig. 3 (retuned LQR for each lambda_h)
	cfg.relaxationRates           = [0.005 0.02 0.05 0.10];
	cfg.runOpenLoop               = true;  % Fig. 1
	cfg.openLoopProfiles          = {'sinusoidal'};   % also 'smooth_step', 'chirp'

	% === Output / checks ===
	cfg.figureVisible  = 'on';             % 'off' for batch runs
	cfg.saveFigures    = true;             % PNG (180 dpi) + FIG in result/figures
	cfg.figureDPI      = 180;
	cfg.selfCheckFlag  = 1;                % mathematical self-checks in P1 (f_self_checks.m)
	cfg.lmiCheckFlag   = 0;                % 1: LMI certificate of P5 (LMI/, needs YALMIP+SDPT3)
	cfg.symbFlag       = symbFlag;
	cfg.problemFolder  = fileparts(mfilename('fullpath'));
	if exist('rootFolder','var')
		cfg.resultFolder = fullfile(rootFolder,'result');
	else
		cfg.resultFolder = fullfile(pwd,'result');
	end

% ===============================================
% 問題設定 / Problem setting (symbolic or analytic, see buildProblem.m)
% ===============================================
	prob = buildProblem(par, cfg);   % f, g, HamF, HamG, Hbar, nu*, canonical field, sigma

	% === Reference-toolbox names (all in the scaled HJ coordinates y) ===
	A = prob.A;                          % A5 = d a_sigma/dy (0)
	B = prob.B;                          % B5 = G_y = e5/S5
	Q = diag(cfg.Qdiag.*cfg.scale.^2);   % Q_y = S*Qz*S   (no factor 1/2)
	R = cfg.Rnu;                         % R_nu
	Qz = diag(cfg.Qdiag);                % physical weight on z - z*

	% === システムの次元(編集しないコト) / System dimension (do not touch) ===
	dim = length(A);

	% === 近似回数・ξの大きさ / Iterations and seed radius (kFlag = xiRFlag = 1) ===
	knum = cfg.manifoldMaxIterations;
	radi = cfg.manifoldRadius;
	tol  = cfg.manifoldIterationRelTol;

	% === 応答性を見るときの初期値と時間 / Initial state and interval of simulation ===
	if isempty(cfg.initialState)
		Initial_State = [0; 0; 0; par.hRelaxEq; cfg.initialRawCommand];
	else
		Initial_State = cfg.initialState(:);
	end
	Time_int = [0, cfg.closedLoopTFinal];
