function [H, stage, J, V, Mu, g, relative] = f_value_colloc(Y, model, des, withJac)
% ==================================================================
% [H, stage, J, V, Mu, g, relative] = f_value_colloc(Y, model, des, withJac)
%
% HJB residual R_H(y) = Hbar(y, grad Vhat(y)) at collocation states Y and its
% Jacobian with respect to the coefficients c (analytic). relative is the
% residual measure rho_H of paper eq. (rho):
%     rho_H = |R_H| / (1e-10 + l + |grad Vhat'*(a_sigma + G_y nu*)|)
% ==================================================================
[V, Mu, Jmu, g] = f_value_eval(Y, model, des, withJac);
[H, stage, term, ~, fs] = f_hamcalcu(Y, Mu, des);
relative = abs(H)./(1e-10 + stage + abs(term));
J = [];
if withJac
    J = zeros(size(Y,2), numel(model.coefficients));
    for j = 1:des.n
        J = J + reshape(Jmu(j,:,:), numel(model.coefficients), [])'.*fs(j,:)';
    end
end
end
