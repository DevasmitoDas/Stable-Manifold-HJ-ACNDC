function [res, J] = f_value_anchor_res(model, manifold, des, cfg)
% ==================================================================
% [res, J] = f_value_anchor_res(model, manifold, des, cfg)
%
% Anchor term of the fit (paper eq. (fit)):
%     sqrt(omega_a) * (grad Vhat(y^m_j) - mu^m_j) / a_j,   a_j = max(1e-7, ||mu^m_j||)
% stacked column by column, and its Jacobian with respect to the coefficients.
% ==================================================================
[~, mu, Jmu] = f_value_eval(manifold.Y, model, des, true);
den = max(1e-7, sqrt(sum(manifold.Mu.^2,1)));
res = reshape((mu - manifold.Mu)./den, [], 1)*sqrt(cfg.manifoldAnchorWeight);
J = [];
if nargout > 1
    nb = numel(model.coefficients);
    for j = 1:des.n
        v = reshape(Jmu(j,:,:), nb, [])./den;
        Jmu(j,:,:) = reshape(v, 1, nb, []);
    end
    J = reshape(permute(Jmu, [1 3 2]), [], nb)*sqrt(cfg.manifoldAnchorWeight);
end
end
