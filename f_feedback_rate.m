function [nu, used, relative, reason] = f_feedback_rate(z, kind, des, model, cfg)
% ==================================================================
% [nu, used, relative, reason] = f_feedback_rate(z, kind, des, model, cfg)
%
% Guarded rate selection at a sample (paper Sec. 4.2):
%  (i)   y = S^{-1}(z - z*) with the RAW u; nu_L = sat_nubar(-K5*y)
%  (ii)  HJ eligible if ||y|| <= 1.6, Vhat and grad Vhat finite, |g| <= 3,
%        rho_H(y) <= 0.05
%  (iii) nu_H = sat_nubar(-G_y'*grad Vhat/(2*Rnu))
%  (iv)  nu_H applied only if grad Vhat'*(a_sigma + G_y nu_H) <= -0.02*y'*Qy*y
%        (eq. (guard)), otherwise nu_L.
% kind = 'LQR' returns nu_L. reason: 0 HJ used / LQR requested, 1 domain,
% 2 non-finite value or |g| > 3, 3 HJB residual, 4 decrease guard.
% ==================================================================
lo = -cfg.uRateMax; hi = cfg.uRateMax;
q = z - des.zStar; y = q./des.scale;
nu = min(max(-des.K*y, lo), hi); used = false; relative = nan; reason = 0;
if ~strcmp(kind,'HJ'), return; end
if norm(y) > cfg.hjDomainRadius, reason = 1; return; end
[V, Mu, ~, g] = f_value_eval(y, model, des, false);
if ~isfinite(V) || abs(g) > cfg.maxLogCorrection || any(~isfinite(Mu)), reason = 2; return; end
[H, stage, term] = f_hamcalcu(y, Mu, des);
relative = abs(H)/(1e-10 + stage + abs(term));
if relative > cfg.hjResidualMax, reason = 3; return; end
lambda = Mu./des.scale; candidate = f_inputcalcu(lambda, des, lo, hi);
qdot = des.scale.*(des.prob.f(y) + des.B*candidate); stateCost = q'*des.Qphysical*q;
if lambda'*qdot > -cfg.hjDecreaseFraction*stateCost + 1e-12, reason = 4; return; end
nu = candidate; used = true;
end
