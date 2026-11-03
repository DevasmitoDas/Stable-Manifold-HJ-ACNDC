function cse = f_simulate_case(kind, par, des, model, cfg, bd, amplitude, Tfinal)
% ==================================================================
% cse = f_simulate_case(kind, par, des, model, cfg, bd, amplitude, Tfinal)
%
% Sampled loop of paper Fig. 2 / Sec. 5:
%   kind = 'LQR' | 'HJ'   : feedback rate from f_feedback_rate at t_k = k*Ts
%   kind = profile name   : open loop, raw rate that tracks sigma^{-1}(profile)
% The rate is held on each interval, u_{k+1} = u_k + Ts*nu_k exactly, and the
% plant is integrated with cfg.rkSubsteps RK4 substeps (sigma(u) evaluated at
% every stage). State cost: RK4-stage quadrature; rate cost: Ts*Rnu*nu_k^2.
%   bd        : disturbance coupling B_d
%   amplitude : scale of the exosystem input w (0: no forcing)
% ==================================================================
pp = par; pp.Bd = bd;
steps = ceil(Tfinal/cfg.Ts); t = (0:steps)'*cfg.Ts; t(end) = Tfinal;
z = zeros(steps+1,5); eta = zeros(steps+1,2); eta(1,:) = cfg.initialEta(:)';
if isempty(cfg.initialState)
    z(1,:) = [0 0 0 pp.hRelaxEq cfg.initialRawCommand];
else
    z(1,:) = cfg.initialState(:)';
end
isFeedback = any(strcmp(kind,{'LQR','HJ'}));
if ~isFeedback
    target0 = f_command_profile(0, kind, pp, Tfinal);
    assert(target0 > cfg.uMin && target0 < cfg.uMax, 'Applied profile starts outside the saturation range.');
    z(1,5) = des.prob.inverseSaturation(target0);
end
nu = zeros(steps,1); used = false(steps,1); hres = nan(steps,1); reasons = zeros(steps,1);
stateCost = zeros(steps+1,1); rateCost = zeros(steps+1,1); status = 'completed'; last = steps+1;
Qz = des.Qphysical;
for k = 1:steps
    interval = t(k+1) - t(k);
    if isFeedback
        [nu(k), used(k), hres(k), reasons(k)] = f_feedback_rate(z(k,:)', kind, des, model, cfg);
    else
        appliedTarget = f_command_profile(t(k+1), kind, pp, Tfinal);
        assert(appliedTarget > cfg.uMin && appliedTarget < cfg.uMax, 'Applied profile is outside the saturation range.');
        target = des.prob.inverseSaturation(appliedTarget);
        nu(k) = min(max((target - z(k,5))/interval, -cfg.uRateMax), cfg.uRateMax);
    end
    state = [z(k,:)'; eta(k,:)']; dt = interval/cfg.rkSubsteps; integral = 0;
    for j = 1:cfg.rkSubsteps
        tj = t(k) + (j-1)*dt;
        k1 = i_state_eq_for_control(tj,      state,        nu(k), pp, des, amplitude);
        s2 = state + dt*k1/2;  k2 = i_state_eq_for_control(tj+dt/2, s2, nu(k), pp, des, amplitude);
        s3 = state + dt*k2/2;  k3 = i_state_eq_for_control(tj+dt/2, s3, nu(k), pp, des, amplitude);
        s4 = state + dt*k3;    k4 = i_state_eq_for_control(tj+dt,   s4, nu(k), pp, des, amplitude);
        q1 = state(1:5) - des.zStar; q2 = s2(1:5) - des.zStar;
        q3 = s3(1:5) - des.zStar;    q4 = s4(1:5) - des.zStar;
        integral = integral + dt/6*(q1'*Qz*q1 + 2*q2'*Qz*q2 + 2*q3'*Qz*q3 + q4'*Qz*q4);
        state = state + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    end
    % Exact held-rate endpoint of the RAW integrator (no amplitude clipping of u)
    state(5) = z(k,5) + interval*nu(k);
    z(k+1,:) = state(1:5)'; eta(k+1,:) = state(6:7)';
    stateCost(k+1) = stateCost(k) + integral; rateCost(k+1) = rateCost(k) + interval*des.R*nu(k)^2;
    if any(~isfinite(state)) || norm((state(1:5) - des.zStar)./des.scale) > cfg.divergenceScaledNorm
        status = 'diverged'; last = k+1; break
    end
end
t = t(1:last); z = z(1:last,:); eta = eta(1:last,:); nu = nu(1:last-1);
used = used(1:last-1); hres = hres(1:last-1); reasons = reasons(1:last-1);
stateCost = stateCost(1:last); rateCost = rateCost(1:last);
residual = eta*pp.Cd'; torque = bd*residual;

cse = struct();
cse.name = kind; cse.controller = kind; cse.status = status;
cse.t = t; cse.z = z; cse.x = z(:,1:4); cse.uRaw = z(:,5);
cse.u = des.prob.saturation(cse.uRaw);                 % APPLIED command sigma(u)
cse.nu = nu;
cse.appliedRateAtSample = des.prob.saturationDerivative(cse.uRaw(1:end-1)).*nu;
cse.eta = eta; cse.d = residual; cse.disturbanceTorque = torque; cse.effectiveAcceleration = torque/pp.J;
cse.hjUsed = used; cse.hjResidual = hres; cse.fallbackReason = reasons;
cse.stateCostIntegral = stateCost; cse.rateCostIntegral = rateCost;
cse.Bd = bd; cse.excitationScale = amplitude; cse.lambda_h = pp.lambda_h; cse.Rnu = des.R;
cse.metrics = f_case_metrics(cse, des, cfg);
end
