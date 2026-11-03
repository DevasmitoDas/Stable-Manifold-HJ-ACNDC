function nu = f_inputcalcu(lambda, des, lo, hi)
% ==================================================================
% nu = f_inputcalcu(lambda, des, lo, hi)
%
% F_INPUTCALCU  optimal raw rate from the PHYSICAL costate lambda (= mu./S)
%     nu* = -G'*lambda/(2*Rnu) = -mu5/(2*Rnu*s_u)          (eq. (nustar))
% evaluated with the symbolically derived stationary rate and, optionally,
% clipped to [lo, hi]. Since Hpre is a convex parabola in nu, the clipped value
% is the exact minimiser of Hpre over [lo, hi] (paper Sec. 4.2 (iii)).
% ==================================================================
if nargin < 3, lo = -inf; hi = inf; end
Mu = des.scale.*lambda;
nu = des.prob.optimalRate(zeros(des.n, size(lambda,2)), Mu);
nu = min(max(nu, lo), hi);
end
