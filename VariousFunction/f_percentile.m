function q = f_percentile(v, p)
% ==================================================================
% q = f_percentile(v, p)
%
% F_PERCENTILE  p-th percentile with linear interpolation (no toolbox needed)
% ==================================================================
v = sort(v(:)); at = 1 + (numel(v)-1)*p/100; lo = floor(at); hi = ceil(at);
q = v(lo) + (at - lo)*(v(hi) - v(lo));
end
