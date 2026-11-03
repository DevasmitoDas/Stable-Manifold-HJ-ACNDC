function src = i_nonlinear_term(z, des)
% ==============================================
% src = i_nonlinear_term(z, des)
%
% Nonlinear remainder [n_s; n_u] of the canonical equations in the
% diagonalised coordinates (s, chi) = T^{-1}(y, mu)  (paper Sec. 3.2):
%     sdot   =  F_c  s  + n_s(s,chi)
%     chidot = -F_c' chi + n_u(s,chi)
% i.e. src = T^{-1}*canonical(T*z) - diag(F_c,-F_c')*z, for many columns z.
% ==============================================
yp = des.T*z;
src = des.Tinv*i_canonical_eq(0, yp, des) - des.diagonalLinear*z;
end
