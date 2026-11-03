function ns = i_highf(z, des)
% ==================================================================
% ns = i_highf(z, des)
%
% Nonlinear term n_s of the stable (x') equation in diagonalised coordinates
% (reference name). z = [s; chi] (2n-by-N).
% ==================================================================
src = i_nonlinear_term(z, des);
ns = src(1:des.n,:);
end
