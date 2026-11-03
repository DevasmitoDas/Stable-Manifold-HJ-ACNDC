function X = i_directforp(nu_src, ker)
% ==================================================================
% X = i_directforp(nu_src, ker)
%
% Backward part of one successive approximation (paper eq. (iter)):
%     chi_{k+1}(t) = -int_t^{T_m} e^{F'(a-t)} n_u(s_k(a), chi_k(a)) da,  chi(T_m) = 0
% evaluated exactly on the mesh ker.t (transposed weights of i_directforx):
%     chi(t_i) = E_i' chi(t_{i+1}) - right_i' n_u(t_i) - left_i' n_u(t_{i+1}).
%   nu_src : n_u along the previous iterate (n-by-nt)
% ==================================================================
nt = numel(ker.t); X = zeros(size(nu_src,1), nt);
for i = nt-1:-1:1
    X(:,i) = ker.E(:,:,i)'*X(:,i+1) - ker.right(:,:,i)'*nu_src(:,i) - ker.left(:,:,i)'*nu_src(:,i+1);
end
end
