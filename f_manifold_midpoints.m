function mid = f_manifold_midpoints(z, source, ker)
% ==================================================================
% mid = f_manifold_midpoints(z, source, ker)
%
% Trajectory at the interval midpoints reconstructed with the half-interval
% weights (forward for s, backward for chi). Used for the "midpoint forcing
% defect" acceptance test of Table 1 (SP_manifold_check.m).
% ==================================================================
nx = size(z,1)/2; nt = size(z,2) - 1; is = 1:nx; iu = nx+(1:nx); mid = zeros(2*nx,nt);
for i = 1:nt
    ns = 0.5*(source(is,i) + source(is,i+1)); un = 0.5*(source(iu,i) + source(iu,i+1));
    mid(is,i) = ker.halfE(:,:,i)*z(is,i) + ker.halfLeft(:,:,i)*source(is,i) + ker.halfRight(:,:,i)*ns;
    mid(iu,i) = ker.halfE(:,:,i)'*z(iu,i+1) - ker.halfRight(:,:,i)'*un - ker.halfLeft(:,:,i)'*source(iu,i+1);
end
end
