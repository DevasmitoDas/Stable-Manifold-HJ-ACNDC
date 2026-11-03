function [V, Mu, JacMu, g] = f_value_eval(Y, model, des, needJacobian)
% ==================================================================
% [V, Mu, JacMu, g] = f_value_eval(Y, model, des, needJacobian)
%
% Positive-definite value ansatz (paper eq. (ansatz), Proposition 2)
%     Vhat(y) = (y'*P5*y) * exp(g(y)),   g(y) = sum_{1<=|alpha|<=3} c_alpha y^alpha
%     grad Vhat = exp(g) * [2*P5*y + (y'*P5*y)*grad g]          ( = mu )
% Vhat(0) = 0, Vhat > 0 for y ~= 0 and Hess Vhat(0) = 2*P5 for EVERY c, and
% the five costate components are the gradient of ONE scalar function.
%
%   Y      : 5-by-np scaled states
%   model  : .exponents (f_poly_exponents), .coefficients (nb-by-1)
%   des    : design (des.P = P5, des.n)
%   JacMu  : d Mu / d c   (n-by-nb-by-np), only if needJacobian
% ==================================================================
if nargin < 4, needJacobian = false; end
[phi, dphi] = f_poly_basis(Y, model.exponents);
g = model.coefficients'*phi; nb = numel(model.coefficients); np = size(Y,2);
dg = zeros(des.n, np);
for j = 1:des.n
    dp = reshape(dphi(j,:,:), nb, np); dg(j,:) = model.coefficients'*dp;
end
V2 = sum(Y.*(des.P*Y),1); base = 2*des.P*Y;
eg = exp(g); V = V2.*eg; Mu = (base + V2.*dg).*eg;
JacMu = [];
if needJacobian
    JacMu = zeros(des.n, nb, np);
    for j = 1:des.n
        dp = reshape(dphi(j,:,:), nb, np);
        Jj = eg.*(phi.*(base(j,:) + V2.*dg(j,:)) + V2.*dp);
        JacMu(j,:,:) = reshape(Jj, 1, nb, np);
    end
end
end
