function model = f_value_fit(des, manifold, Ycoll, cfg)
% ==================================================================
% model = f_value_fit(des, manifold, Ycoll, cfg)
%
% Fit of the positive-definite value ansatz (paper eq. (fit)):
%   min_c  sum_i [w_i R_H(y_i)]^2 + omega_a sum_j ||grad Vhat(y^m_j) - mu^m_j||^2/a_j^2
%          + lambda_c ||c||^2
%   w_i = 1/(0.03 + l_i) frozen at each step, Gauss-Newton with analytic Jacobians,
%   backtracking line search that also enforces |g| <= 3 on the collocation set.
%
%   manifold : .Y, .Mu  manifold anchors (scaled state / scaled costate)
%   Ycoll    : collocation states (built in SP_value_fit.m)
% ==================================================================
model = struct('exponents', f_poly_exponents(cfg.valueCorrectionDegree, des.n));
model.coefficients = zeros(size(model.exponents,1), 1);
model.fitHistory = [];
nb = numel(model.coefficients);
for it = 1:cfg.fitMaxIterations
    [H, stage, J] = f_value_colloc(Ycoll, model, des, true);
    weight = 1./(cfg.fitWeightFloor + stage(:));
    r = weight.*H(:);
    Aw = weight.*J;
    [ra, Ja] = f_value_anchor_res(model, manifold, des, cfg);
    delta = -[Aw; Ja; sqrt(cfg.fitRidge)*eye(nb)] \ ...
        [r; ra; sqrt(cfg.fitRidge)*model.coefficients];
    objective = sum(r.^2) + sum(ra.^2) + cfg.fitRidge*sum(model.coefficients.^2);
    accepted = false; step = 1;
    for ls = 1:12
        trial = model; trial.coefficients = model.coefficients + step*delta;
        [Ht, ~, ~, ~, ~, gt] = f_value_colloc(Ycoll, trial, des, false);
        if all(isfinite(Ht)) && max(abs(gt)) <= cfg.maxLogCorrection
            rat = f_value_anchor_res(trial, manifold, des, cfg);
            ot = sum((weight.*Ht(:)).^2) + sum(rat.^2) + cfg.fitRidge*sum(trial.coefficients.^2);
            if ot < objective
                model.coefficients = trial.coefficients; accepted = true; break
            end
        end
        step = step/2;
    end
    model.fitHistory(end+1,:) = [it objective accepted step norm(delta)];
    fprintf('  HJB fit %d: weighted objective=%g, step accepted=%d\n', it, objective, accepted);
    if ~accepted || norm(step*delta) < 1e-6
        break
    end
end
model.trainingPoints = size(Ycoll,2);
model.stateDimension = des.n;
model.coefficientCount = nb;
model.description = ['Vhat=(y''*P5*y)*exp(polynomial g(y)); exact quadratic Hessian; ' ...
    'local manifold anchors + approximate HJB collocation.'];
end
