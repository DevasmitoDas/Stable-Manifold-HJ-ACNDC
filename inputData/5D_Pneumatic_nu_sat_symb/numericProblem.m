function prob = numericProblem(par, cfg)
%NUMERICPROBLEM Analytic (toolbox-free) form of the problem defined in symbSetting.m.
%   prob = numericProblem(par, cfg)   (called by buildProblem.m when symbFlag = 0
%                                       or when the Symbolic Math Toolbox is missing)
%
% Scaled error y = S^{-1}(z - z*), q = S*y, scaled costate mu = S*lambda.
%   sigma(u)  = u_m + u_h*tanh((u-u_m)/u_h),  sigma'(u) = 1 - tanh^2(.)
%   r         = (-q3 + ku*(sigma(u*+q5) - sigma(u*)))/taup              (eq. (press))
%   hdot      = (Ah - gamma*s_eh(h))*r - beta*h*s_er(r) - lambda_h*q4    (eq. (hyst))
%   F(q)      = [q2; (-k q1 - c q2 + bp q3 + kh q4)/J; r; hdot; 0],  f(y) = S^{-1}F(Sy)
%   Hbar      = mu'*f(y) + y'*Qy*y - (mu5/S5)^2/(4*Rnu)                  (eq. (Hbar))
%   nu*       = -mu5/(2*Rnu*S5)                                          (eq. (nustar))
%   lambdadot = -(dF/dq)'*lambda - 2*Qz*q, in particular
%   lambda5dot = -(ku/taup)*sigma'(u)*(lambda3 + Dr*lambda4) - 2*q_u*(u-u*)  (eq. (lam5))
% These are exactly the expressions that symbSetting.m derives symbolically; the
% two routes agree to about 1e-13 (relative) on f, Hbar and the canonical field.

m = model_constants(par, cfg);

prob = struct();
prob.n = 5;
prob.nInput = 1;
prob.zStar = m.zStar;
prob.uAppliedStar = m.uaStar;

prob.saturation           = @(u)  f_SatTanh(u, cfg.uMin, cfg.uMax);
prob.saturationDerivative = @(u)  f_DSatTanh(u, cfg.uMin, cfg.uMax);
prob.inverseSaturation    = @(ua) f_SatTanh_inv(ua, cfg.uMin, cfg.uMax);

Bscaled = m.B;
prob.f  = @(Y) drift(Y, m);
prob.g  = @(Y) Bscaled;
prob.driftJacobian  = @(y) drift_jacobian(y, m);
prob.optimalRate    = @(Y, Mu) -Mu(5,:)/(2*m.R*m.S(5));
prob.preHamiltonian = @(Y, Mu, nu) sum(Mu.*(drift(Y,m) + Bscaled*nu),1) ...
                                   + sum(Y.*(m.Qy*Y),1) + m.R*nu.^2;
prob.dPreDnu        = @(Y, Mu, nu) Mu(5,:)/m.S(5) + 2*m.R*nu;
prob.hamiltonian    = @(Y, Mu) sum(Mu.*drift(Y,m),1) + sum(Y.*(m.Qy*Y),1) ...
                                   - Mu(5,:).^2/(4*m.R*m.S(5)^2);
prob.canonical      = @(Y, Mu) canonical_field(Y, Mu, m);

Alin = drift_jacobian(zeros(5,1), m);
Dlin = (Bscaled*Bscaled')/(2*m.R);
Hlin = [Alin, -Dlin; -2*m.Qy, -Alin'];
prob.A = Alin;
prob.B = Bscaled;
prob.canonicalLinear = Hlin;
prob.HamF = @(Y, Mu) nonlinear_part(Y, Mu, m, Hlin, 1:5);
prob.HamG = @(Y, Mu) nonlinear_part(Y, Mu, m, Hlin, 6:10);

prob.inputCurvature = 2*m.R;
prob.sourceFolder = fileparts(mfilename('fullpath'));
prob.status = 'Analytic saturation-aware formulas (numericProblem.m) active';
prob.generatedFolder = '';
prob.generated = false;
prob.stationarityIdentityVerified = true;   % nu* is the exact root of dPreDnu
end

% ------------------------------------------------------------------------------
function m = model_constants(par, cfg)
m = struct('J',par.J,'k',par.k,'c',par.c,'bp',par.bp,'kh',par.kh, ...
    'taup',par.taup,'ku',par.ku,'Ah',par.Ah,'beta',par.beta,'gamma',par.gamma, ...
    'lambda_h',par.lambda_h,'hr',par.hRelaxEq,'eh',par.eps_abs_h,'er',par.eps_abs_rate, ...
    'uMin',cfg.uMin,'uMax',cfg.uMax);
m.S  = cfg.scale(:);
m.R  = cfg.Rnu;
m.Qz = diag(cfg.Qdiag(:));
m.Qy = diag(cfg.Qdiag(:).*m.S.^2);
pStar = (par.k*cfg.thetaStar - par.kh*par.hRelaxEq)/par.bp;
m.uaStar = pStar/par.ku;
if ~(m.uaStar > cfg.uMin && m.uaStar < cfg.uMax)
    error('HJ:Equilibrium','The target APPLIED command must lie strictly inside the saturation limits.');
end
m.uStar = f_SatTanh_inv(m.uaStar, cfg.uMin, cfg.uMax);   % raw u* = sigma^{-1}(0.6) = 0.60137
m.zStar = [cfg.thetaStar; 0; pStar; par.hRelaxEq; m.uStar];
m.B = [0; 0; 0; 0; 1]./m.S;                              % G_y = S^{-1} e5
end

function [r, hd, Dr, Dh, ds] = rate_terms(q, m)
% Pressure rate r, hysteresis rate hd and the partial derivatives
% Dr = d(hd)/dr, Dh = d(hd)/dh, ds = sigma'(u).
uraw = m.uStar + q(5,:);
ua = f_SatTanh(uraw, m.uMin, m.uMax);
ds = f_DSatTanh(uraw, m.uMin, m.uMax);
r  = (-q(3,:) + m.ku*(ua - m.uaStar))/m.taup;
h  = m.hr + q(4,:);
sh = f_smooth_abs(h, m.eh);
sr = f_smooth_abs(r, m.er);
a  = m.Ah - m.gamma*sh;
hd = a.*r - m.beta*h.*sr - m.lambda_h*q(4,:);
Dr = a - m.beta*h.*r./hypot(r, m.er);
Dh = -m.gamma*h./hypot(h, m.eh).*r - m.beta*sr - m.lambda_h;
end

function Fy = drift(Y, m)
q = m.S.*Y;
[r, hd] = rate_terms(q, m);
Fq = [q(2,:); ...
      (-m.k*q(1,:) - m.c*q(2,:) + m.bp*q(3,:) + m.kh*q(4,:))/m.J; ...
      r; ...
      hd; ...
      zeros(1,size(q,2))];
Fy = Fq./m.S;
end

function Ay = drift_jacobian(y, m)
q = m.S.*y(:);
[~, ~, Dr, Dh, ds] = rate_terms(q, m);
Jq = zeros(5);
Jq(1,2) = 1;
Jq(2,1:4) = [-m.k, -m.c, m.bp, m.kh]/m.J;
Jq(3,3) = -1/m.taup;
Jq(3,5) = m.ku*ds/m.taup;
Jq(4,3) = -Dr/m.taup;
Jq(4,4) = Dh;
Jq(4,5) = Dr*m.ku*ds/m.taup;
Ay = (Jq.*m.S')./m.S;          % S^{-1} * Jq * S
end

function out = canonical_field(Y, Mu, m)
% [ydot; mudot] = [dHbar/dmu; -dHbar/dy]   (paper eq. (canon))
q   = m.S.*Y;
lam = Mu./m.S;
[r, hd, Dr, Dh, ds] = rate_terms(q, m);
nu  = -Mu(5,:)/(2*m.R*m.S(5));
zdot = [q(2,:); ...
        (-m.k*q(1,:) - m.c*q(2,:) + m.bp*q(3,:) + m.kh*q(4,:))/m.J; ...
        r; ...
        hd; ...
        nu];
ldot = [ m.k/m.J*lam(2,:); ...
        -lam(1,:) + m.c/m.J*lam(2,:); ...
        -m.bp/m.J*lam(2,:) + lam(3,:)/m.taup + Dr.*lam(4,:)/m.taup; ...
        -m.kh/m.J*lam(2,:) - Dh.*lam(4,:); ...
        -m.ku/m.taup*ds.*(lam(3,:) + Dr.*lam(4,:))] - 2*m.Qz*q;
out = [zdot./m.S; m.S.*ldot];
end

function out = nonlinear_part(Y, Mu, m, Hlin, rows)
% HamF / HamG: canonical field minus its linearisation Ham*[y;mu]
full = canonical_field(Y, Mu, m) - Hlin*[Y; Mu];
out = full(rows,:);
end
