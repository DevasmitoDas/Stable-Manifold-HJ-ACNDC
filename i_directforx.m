function S = i_directforx(xi, ns, ker)
% ==================================================================
% S = i_directforx(xi, ns, ker)
%
% Forward part of one successive approximation (paper eq. (iter)):
%     s_{k+1}(t) = e^{F t} xi + int_0^t e^{F(t-a)} n_s(s_k(a), chi_k(a)) da
% evaluated exactly on the mesh ker.t for a forcing that is linear on each
% interval:  s(t_{i+1}) = E_i s(t_i) + left_i n_s(t_i) + right_i n_s(t_{i+1}).
%   xi : seed (n-by-1),  ns : n_s along the previous iterate (n-by-nt)
% (In the reference toolbox this was the ode45 integrand; here the integral
%  is computed with exact exponential weights, f_manifold_kernels.m.)
% ==================================================================
nt = numel(ker.t); S = zeros(numel(xi), nt); S(:,1) = xi;
for i = 1:nt-1
    S(:,i+1) = ker.E(:,:,i)*S(:,i) + ker.left(:,:,i)*ns(:,i) + ker.right(:,:,i)*ns(:,i+1);
end
end
