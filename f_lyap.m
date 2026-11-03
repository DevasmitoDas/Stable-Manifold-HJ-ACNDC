function L = f_lyap(F, D)
% ==================================================================
% L = f_lyap(F, D)
%
% F_LYAP  symmetric solution of  F*L + L*F' = D  (no Control System Toolbox)
%   Paper Sec. 3.2: F_c*L + L*F_c' = D, D = B5*B5'/(2*Rnu).
%   Toolbox relation: L = lyap(F, -Rbar).
% ==================================================================
nx = size(F,1);
M = kron(eye(nx),F) + kron(F,eye(nx));
L = reshape(M\D(:), nx, nx);
L = (L + L')/2;
end
