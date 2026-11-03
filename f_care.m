function P = f_care(A, B, Q, R)
% ==================================================================
% P = f_care(A, B, Q, R)
%
% F_CARE  stabilising solution of the CARE (no Control System Toolbox)
%     A'*P + P*A - P*B*R^{-1}*B'*P + Q = 0      (paper: P5, V2 = y'*P5*y)
% from the stable invariant subspace of [A -B*B'/R; -Q -A'].
% Toolbox relation: care(A,B,2*Q,2*R) = 2*P5 = Gamma (costate gain p = P x).
% ==================================================================
n = size(A,1); H = [A -(B*B')/R; -Q -A'];
[U,E] = eig(H); ix = find(real(diag(E)) < -1e-9);
if numel(ix) ~= n
    error('HJ:Riccati','Hamiltonian has no separated n-dimensional stable subspace.');
end
U = U(:,ix);
if rcond(U(1:n,:)) < 1e-12
    error('HJ:Riccati','Ill-conditioned Riccati stable subspace.');
end
Pc = U(n+1:end,:)/U(1:n,:);
if norm(imag(Pc),'fro') > 1e-7*max(1,norm(Pc,'fro'))
    error('HJ:Riccati','Riccati solution is not real.');
end
P = real(Pc); P = (P + P')/2;
res = A'*P + P*A - P*(B*B')*P/R + Q;
assert(norm(res,'fro') < 1e-6*max(1,norm(Q,'fro')), 'CARE residual too large.');
assert(min(eig(P)) > 0, 'Riccati matrix must be positive definite.');
end
