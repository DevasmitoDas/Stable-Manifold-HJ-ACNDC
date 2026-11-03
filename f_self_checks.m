function checks = f_self_checks(par, des, cfg)
% ==================================================================
% checks = f_self_checks(par, des, cfg)
%
% Mathematical self-checks run in P1 (paper Sec. 5.1 quotes several of them):
%  - symplecticity of T and block diagonalisation T^{-1}*Ham*T = diag(F,-F')
%  - linearisation of the augmented plant (finite differences)
%  - stationarity of nu*, bounded minimiser over [-nubar, nubar], affineness in nu
%  - symbolic substitution  Hpre(y,mu,nu*) = Hbar(y,mu)
%  - chain rule: lambda5dot contains sigma'(u)   (eq. (lam5)), far from u*
%  - saturation bounds 0 < sigma < 1, 0 < sigma' <= 1, sigma(u*) = 0.6
%  - canonical field = symplectic gradient of Hbar (finite differences)
%  - exact exponential weights of the integral iteration
%  - gradient of the value ansatz and the collocation Jacobian (finite differences)
% Any failure stops the run (assert), as in the Revision-5 code.
% ==================================================================
nx = des.n; d = des; c = cfg;
checks = struct('equilibriumResidual',d.equilibriumResidual,'controllabilityRank',d.controllabilityRank, ...
    'riccatiResidual',d.riccatiResidual,'stateDimension',nx,'canonicalDimension',2*nx);

Jsym = [zeros(nx) eye(nx); -eye(nx) zeros(nx)];
checks.symplecticResidual = norm(d.T'*Jsym*d.T - Jsym,'fro');
checks.symplecticDiagonalizationResidual = norm(d.Tinv*d.canonicalLinear*d.T - d.diagonalLinear,'fro');
assert(checks.symplecticResidual < 1e-6 && checks.symplecticDiagonalizationResidual < 1e-6, ...
    'Canonical transform is inconsistent.');

% --- linearisation of the augmented plant ---
hh = 1e-7; An = zeros(nx);
Bn = (shifted_rhs(zeros(nx,1), hh, d) - shifted_rhs(zeros(nx,1), -hh, d))/(2*hh);
for j = 1:nx
    dq = zeros(nx,1); dq(j) = hh;
    An(:,j) = (shifted_rhs(dq, 0, d) - shifted_rhs(-dq, 0, d))/(2*hh);
end
checks.linearizationError = max(norm(An - d.Aphysical,'fro'), norm(Bn - d.Bphysical));
assert(checks.linearizationError < 1e-5, 'Augmented linearisation disagrees with the plant.');

% --- rate minimisation (Proposition 1(a)) ---
Y  = 0.2*(2*f_halton(20,[2 3 5 7 11],101) - 1);
Mu = 2*f_halton(20,[13 17 19 23 29],151) - 1;
worst = 0; stationarity = 0;
for i = 1:size(Y,2)
    q = d.scale.*Y(:,i); lambda = Mu(:,i)./d.scale;
    nr = f_inputcalcu(lambda, d);
    stationarity = max(stationarity, abs(lambda(5) + 2*d.R*nr));
    u = c.uMin + (c.uMax - c.uMin)*(i-1)/(size(Y,2)-1);
    lo = -c.uRateMax; hi = c.uRateMax; nb = f_inputcalcu(lambda, d, lo, hi);
    nuGrid = linspace(lo, hi, 1001); values = lambda(5)*nuGrid + d.R*nuGrid.^2;
    worst = max(worst, lambda(5)*nb + d.R*nb^2 - min(values));
    applied = d.prob.saturation(u + c.Ts*nb);
    assert(applied >= c.uMin - 1e-12 && applied <= c.uMax + 1e-12, 'Applied saturation bound failed.');
    f0 = shifted_rhs(q, 0, d); f1 = shifted_rhs(q, 1.7, d);
    assert(norm(f1 - f0 - 1.7*d.G) < 1e-12, 'Model is not affine in nu.');
end

% --- symbolic preH -> derivative -> stationary rate -> substitution ---
nuSymbolic = d.prob.optimalRate(Y, Mu);
preH = d.prob.preHamiltonian(Y, Mu, nuSymbolic);
Hbar = d.prob.hamiltonian(Y, Mu);
checks.preHamiltonianSubstitutionError = max(abs(preH - Hbar));
checks.symbolicStationarityError = max(abs(d.prob.dPreDnu(Y, Mu, nuSymbolic)));
assert(checks.preHamiltonianSubstitutionError < 1e-10 && checks.symbolicStationarityError < 1e-10, ...
    'Symbolic derivative/substitution identity failed.');

% --- sigma'(u) in the FIFTH costate equation, well outside the linear region ---
chainError = 0;
for raw = [-2 -0.4 0.5 1.4 3]
    q = [.02; -.01; .03; .04; raw - d.uRawStar];
    yy = q./d.scale; mm = [.03; .04; -.02; .05; .06]; ll = mm./d.scale;
    vec = d.prob.canonical(yy, mm); lambdaDot = vec(nx+1:end)./d.scale;
    rr = shifted_rhs(q, 0, d); rr = rr(3);
    Dr = hysteresis_Dr(d.zStar(4) + q(4), rr, par);
    expected = -par.ku/par.taup*d.prob.saturationDerivative(raw)*(ll(3) + Dr*ll(4)) ...
        - 2*d.Qphysical(5,5)*q(5);
    chainError = max(chainError, abs(lambdaDot(5) - expected));
end
checks.saturationCostateChainRuleError = chainError;
assert(chainError < 1e-9, 'The costate equation lost the saturation derivative.');

% --- saturation (Proposition 1(c)) ---
raw = linspace(-10, 10, 401); applied = d.prob.saturation(raw);
slope = d.prob.saturationDerivative(raw);
checks.saturationBoundsPassed = all(applied >= c.uMin - 1e-12 & applied <= c.uMax + 1e-12);
assert(checks.saturationBoundsPassed && all(slope >= -1e-12 & slope <= 1 + 1e-12), ...
    'Smooth saturation bounds/slope are inconsistent.');
assert(abs(d.prob.saturation(d.uRawStar) - d.uAppliedStar) < 1e-12, ...
    'Inverse-saturation equilibrium is inconsistent.');
checks.saturationSlopeAtEquilibrium = d.prob.saturationDerivative(d.uRawStar);   % 0.96
checks.inputStationarityError = stationarity; checks.boundedMinimizerExcess = worst;
checks.inputCurvature = 2*d.R;
assert(stationarity < 1e-12 && worst < 1e-10, 'Quadratic rate minimisation failed.');

% --- canonical field ---
checks.canonicalZeroResidual = norm(i_nonlinear_term(zeros(2*nx,1), d));
assert(checks.canonicalZeroResidual < 1e-12, 'Canonical origin is not an equilibrium.');
y = [0.07; -0.03; 0.02; 0.01; 0.04]; mu = [0.1; 0.02; -0.01; -0.4; 0.03]; yp = [y; mu];
actual = i_canonical_eq(0, yp, d); grad = zeros(2*nx,1); h = 2e-7;
for j = 1:2*nx
    v = zeros(2*nx,1); v(j) = h;
    grad(j) = (point_hamiltonian(yp + v, d) - point_hamiltonian(yp - v, d))/(2*h);
end
expected = [grad(nx+1:end); -grad(1:nx)];
checks.canonicalGradientError = norm(actual - expected)/max(1, norm(expected));
assert(checks.canonicalGradientError < 1e-5, 'Five-state canonical Hamiltonian derivatives disagree.');

% --- exact exponential weights ---
dt = 0.01; [E, left, right] = f_semigroup_weights(d.F, dt);
n0 = linspace(-.2, .4, nx)'; n1 = linspace(.5, -.3, nx)'; tg = linspace(0, dt, 1001);
fw = zeros(nx, numel(tg)); bw = fw;
for k = 1:numel(tg)
    a = tg(k); na = (1 - a/dt)*n0 + (a/dt)*n1;
    fw(:,k) = expm(d.F*(dt - a))*na; bw(:,k) = expm(d.F'*a)*na;
end
checks.integralKernelError = max(norm(left*n0 + right*n1 - trapz(tg, fw, 2)), ...
    norm(right'*n0 + left'*n1 - trapz(tg, bw, 2)));
assert(checks.integralKernelError < 1e-8 && max(abs(eig(E))) < 1, 'Forward/backward semigroup error.');

% --- value ansatz derivatives ---
vm = struct('exponents', f_poly_exponents(3, nx));
vm.coefficients = .01*sin((1:size(vm.exponents,1))');
[~, vg] = f_value_eval(y, vm, d, false); gn = zeros(nx,1);
for j = 1:nx
    v = zeros(nx,1); v(j) = h;
    gn(j) = (f_value_eval(y + v, vm, d, false) - f_value_eval(y - v, vm, d, false))/(2*h);
end
checks.cubicCoefficientCount = numel(vm.coefficients);
checks.valueGradientError = norm(vg - gn)/max(1, norm(vg));
assert(checks.cubicCoefficientCount == 55 && checks.valueGradientError < 1e-6, ...
    'Five-variable scalar fit derivative failed.');
[~, ~, Aj] = f_value_colloc(y, vm, d, true); Jn = zeros(1, numel(vm.coefficients));
for k = 1:numel(vm.coefficients)
    vp = vm; vn = vm;
    vp.coefficients(k) = vp.coefficients(k) + h; vn.coefficients(k) = vn.coefficients(k) - h;
    Jn(k) = (f_value_colloc(y, vp, d, false) - f_value_colloc(y, vn, d, false))/(2*h);
end
checks.collocationJacobianError = norm(Aj - Jn)/max(1, norm(Jn));
assert(checks.collocationJacobianError < 1e-5, 'HJB parameter Jacobian failed its finite-difference test.');

% --- cold start (paper Sec. 5) ---
if isempty(c.initialState)
    initial = [0; 0; 0; par.hRelaxEq; c.initialRawCommand];
else
    initial = c.initialState(:);
end
checks.coldStartScaledRadius = norm((initial - d.zStar)./d.scale);     % 1.733 > 1.6
checks.initialAppliedCommand = d.prob.saturation(initial(5));          % sigma(0) = 0.1192
fprintf(['Self-checks passed: saturation, preH substitution, sigma''(u) chain rule, CARE, ' ...
    'canonical gradients, kernels and value derivatives.\n']);
end

% ------------------------------------------------------------------------------
function f = shifted_rhs(q, nu, d)
% augmented drift in physical error coordinates: qdot = S*(a_sigma(q/S) + G_y*nu)
f = d.scale.*(d.prob.f(q./d.scale) + d.B*nu);
end

function Dr = hysteresis_Dr(h, r, p)
% Dr = d(hdot)/dr = Ah - gamma*s_eh(h) - beta*h*r/sqrt(r^2+eps_r^2)
sh = f_smooth_abs(h, p.eps_abs_h);
Dr = p.Ah - p.gamma*sh - p.beta*h.*r./hypot(r, p.eps_abs_rate);
end

function H = point_hamiltonian(yp, d)
H = f_hamcalcu(yp(1:d.n), yp(d.n+1:2*d.n), d);
end
