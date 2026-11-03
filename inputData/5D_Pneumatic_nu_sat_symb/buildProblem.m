function prob = buildProblem(par, cfg)
%BUILDPROBLEM Problem functions of the saturation-aware pneumatic HJ problem.
%   prob = buildProblem(par, cfg)
%
%   cfg.symbFlag = 1 : symbSetting.m -> symbCompile.m (Symbolic Math Toolbox).
%                      This is the route used for the paper: the pre-Hamiltonian
%                      is formed, d/dnu = 0 is solved, nu* is substituted and the
%                      canonical field is obtained by differentiation.
%   cfg.symbFlag = 0 : numericProblem.m, the same formulas written out by hand
%                      (no toolbox). Used automatically if the toolbox is missing.
%
%   Both routes return the same fields (vectorised over columns, y and mu are 5-by-N):
%     n, nInput, zStar, uAppliedStar, A, B, canonicalLinear,
%     f(y), g(y), HamF(y,mu), HamG(y,mu), hamiltonian(y,mu), optimalRate(y,mu),
%     preHamiltonian(y,mu,nu), dPreDnu(y,mu,nu), canonical(y,mu), driftJacobian(y),
%     saturation(u), saturationDerivative(u), inverseSaturation(ua),
%     inputCurvature, route, status, generatedFolder, generated,
%     stationarityIdentityVerified

useSymbolic = false;
if isfield(cfg,'symbFlag') && cfg.symbFlag == 1
    useSymbolic = symbolic_toolbox_available();
    if ~useSymbolic
        warning('HJ:SymbolicFallback', ['Symbolic Math Toolbox is not available: ' ...
            'using the equivalent analytic formulas of numericProblem.m.']);
    end
end

if useSymbolic
    prob = symbCompile(par, cfg);
    prob.route = 'symbolic';
else
    prob = numericProblem(par, cfg);
    prob.route = 'analytic';
end
end

function ok = symbolic_toolbox_available()
ok = false;
if exist('OCTAVE_VERSION','builtin')
    return
end
try
    ok = (exist('sym','file') == 2) && license('test','Symbolic_Toolbox');
catch
    ok = false;
end
end
