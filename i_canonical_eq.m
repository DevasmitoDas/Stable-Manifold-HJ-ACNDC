function dxp = i_canonical_eq(t, xp, des) %#ok<INUSL>
% ==============================
% DXP = I_CANONICAL_EQ(t, xp, des)
%
% Hamilton canonical equations of the minimised Hamiltonian (paper eq. (canon))
%     d/dt [y; mu] = [dHbar/dmu; -dHbar/dy] = Ham*[y; mu] + [HamF(y,mu); HamG(y,mu)]
% with Ham = [A5 -D; -2Qy -A5'] (des.canonicalLinear). The fifth costate
% equation contains sigma'(u) (eq. (lam5)). xp = [y; mu] may hold many columns.
% The field is the one derived in symbSetting.m (or numericProblem.m).
% ==============================
n = des.n;
dxp = des.prob.canonical(xp(1:n,:), xp(n+1:2*n,:));
end
