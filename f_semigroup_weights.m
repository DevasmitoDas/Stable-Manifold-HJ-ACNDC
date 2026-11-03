function [E, left, right] = f_semigroup_weights(F, dt)
% ==================================================================
% [E, left, right] = f_semigroup_weights(F, dt)
%
% Exact exponential weights for a forcing n(t) interpolated linearly from n0
% to n1 on [0,dt]:
%     z(dt) = E*z(0) + left*n0 + right*n1,   E = expm(F*dt)
% obtained from one augmented matrix exponential (no inverse of F needed).
% ==================================================================
n = size(F,1); zero = zeros(n);
M = expm([F eye(n) zero; zero zero eye(n); zero zero zero]*dt);
E = M(1:n,1:n); total = M(1:n,n+1:2*n); right = M(1:n,2*n+1:3*n)/dt;
left = total - right;
end
