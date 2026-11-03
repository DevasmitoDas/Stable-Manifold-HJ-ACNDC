function nu_term = i_highg(z, des)
% ==================================================================
% nu_term = i_highg(z, des)
%
% Nonlinear term n_u of the unstable (p') equation in diagonalised
% coordinates (reference name). z = [s; chi] (2n-by-N).
% ==================================================================
src = i_nonlinear_term(z, des);
nu_term = src(des.n+1:2*des.n,:);
end
