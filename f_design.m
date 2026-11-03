function d = f_design(prob, cfg)
% ==================================================================
% d = f_design(prob, cfg)
%
% Local (k = 0) design and block diagonalisation (paper Sec. 3.2):
%   A5, B5        linearisation of the scaled augmented drift
%   P5            CARE A5'P5 + P5A5 - P5B5 Rnu^{-1} B5'P5 + Qy = 0
%   K5            = Rnu^{-1} B5'P5,  Gamma = 2*P5,  D = B5B5'/(2Rnu)
%   F_c           = A5 - B5*K5 = A5 - D*Gamma  (Hurwitz)
%   L             F_c*L + L*F_c' = D
%   T             = [I L; Gamma Gamma*L+I],  T^{-1}*Ham*T = diag(F_c, -F_c')
%   T_m           = max(40, 8/|max Re lambda(F_c)|)   (manifold horizon)
% The fields keep the names of the Revision-5 code that produced the paper.
% ==================================================================
nx = prob.n;
d.n = prob.n; d.nInput = prob.nInput; d.prob = prob;
d.scale = cfg.scale(:); d.Qphysical = diag(cfg.Qdiag(:)); d.R = cfg.Rnu; d.Rnu = cfg.Rnu;
d.zStar = prob.zStar(:); d.xStar = d.zStar(1:4); d.uStar = d.zStar(5);
d.uRawStar = d.uStar; d.uAppliedStar = prob.uAppliedStar; d.nuStar = 0;
d.A = prob.A; d.B = prob.B(:); d.Gy = d.B;
d.Aphysical = (d.A.*d.scale)./d.scale'; d.Bphysical = d.B.*d.scale;
d.G = d.Bphysical; d.Q = diag(cfg.Qdiag(:).*d.scale.^2);
d.P = f_care(d.A, d.B, d.Q, d.R);
d.Gamma = 2*d.P; d.K = (d.B'*d.P)/d.R; d.Kphysical = d.K./d.scale';
d.D = (d.B*d.B')/(2*d.R); d.F = d.A - d.D*d.Gamma; d.poles = eig(d.F);
assert(all(real(d.poles) < 0), 'Saturation-aware local Riccati matrix is not Hurwitz.');
d.L = f_lyap(d.F, d.D);
d.T = [eye(nx) d.L; d.Gamma d.Gamma*d.L+eye(nx)]; d.Tinv = d.T\eye(2*nx);
d.canonicalLinear = prob.canonicalLinear;
expected = [d.A -d.D; -2*d.Q -d.A'];
assert(norm(d.canonicalLinear - expected,'fro') < 1e-8, ...
    'Minimised Hamiltonian linearisation disagrees with the CARE convention.');
d.diagonalLinear = blkdiag(d.F, -d.F');
d.horizon = max(cfg.manifoldMinimumHorizon, cfg.manifoldTailExponents/(-max(real(d.poles))));
d.equilibriumResidual = norm(d.scale.*prob.f(zeros(nx,1)));
W = d.B; for k = 1:nx-1, W = [W d.A^k*d.B]; end %#ok<AGROW>
d.controllabilityRank = rank(W);
d.riccatiResidual = norm(d.A'*d.P + d.P*d.A - d.P*(d.B*d.B')*d.P/d.R + d.Q,'fro');
d.diagonalizationResidual = norm(d.Tinv*d.canonicalLinear*d.T - d.diagonalLinear,'fro');
assert(d.equilibriumResidual < 1e-11 && d.controllabilityRank == nx, ...
    'Invalid saturated-command equilibrium or local controllability.');
end
