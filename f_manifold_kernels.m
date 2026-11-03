function ker = f_manifold_kernels(d, c, refinement)
% ==================================================================
% ker = f_manifold_kernels(des, cfg, refinement)
%
% Non-uniform mesh of [0, T_m] (paper Sec. 3.2): steps 0.01 s on [0,4],
% 0.1 s on [4,30] and 0.5 s up to T_m (divided by 2^refinement), and the cached
% exact semigroup weights on every interval and half interval.
%   ker.t                    mesh (row vector)
%   ker.E/left/right(:,:,i)  weights on [t_i, t_{i+1}]
%   ker.halfE/halfLeft/halfRight(:,:,i)  weights on the half interval
%   ker.linear(:,:,i)        expm(F*t_i)  (k = 0 solution s_0(t) = e^{Ft} xi)
% ==================================================================
factor = 2^refinement; last = d.horizon; first = min(4,last); second = min(30,last);
t = unique([0:c.manifoldFastStep/factor:first, first, ...
    first:c.manifoldMiddleStep/factor:second, second, ...
    second:c.manifoldSlowStep/factor:last, last]);
n = numel(t);
ker = struct('t',t,'E',zeros(d.n,d.n,n-1),'left',zeros(d.n,d.n,n-1), ...
    'right',zeros(d.n,d.n,n-1),'halfE',zeros(d.n,d.n,n-1), ...
    'halfLeft',zeros(d.n,d.n,n-1),'halfRight',zeros(d.n,d.n,n-1),'linear',zeros(d.n,d.n,n));
ker.linear(:,:,1) = eye(d.n); cachedDt = []; cached = cell(1,0);
for i = 1:n-1
    dt = t(i+1) - t(i); hit = find(abs(cachedDt - dt) < 1e-11, 1);
    if isempty(hit)
        [E,left,right] = f_semigroup_weights(d.F, dt);
        [Eh,lh,rh] = f_semigroup_weights(d.F, dt/2);
        cachedDt(end+1) = dt; cached{end+1} = {E,left,right,Eh,lh,rh}; hit = numel(cachedDt); %#ok<AGROW>
    end
    a = cached{hit}; ker.E(:,:,i) = a{1}; ker.left(:,:,i) = a{2}; ker.right(:,:,i) = a{3};
    ker.halfE(:,:,i) = a{4}; ker.halfLeft(:,:,i) = a{5}; ker.halfRight(:,:,i) = a{6};
    ker.linear(:,:,i+1) = ker.E(:,:,i)*ker.linear(:,:,i);
end
end
