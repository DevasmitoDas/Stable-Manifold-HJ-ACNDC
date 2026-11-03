function f_write_metrics_csv(file, cases)
% ==================================================================
% f_write_metrics_csv(file, cases)
%
% One CSV row per simulated case (cell array of f_simulate_case outputs).
% ==================================================================
names = {'Case','Controller','Status','Bd','LambdaH','Rnu','ThetaRMSE','TailThetaRMSE', ...
    'FinalAugmentedStateNorm','FinalPhysicalStateNorm','FinalAppliedCommandError', ...
    'FinalRawCommandError','FinalPressureError','FinalHysteresisError','SettlingTime', ...
    'StateCost','RateCost','FiniteHorizonCost','MinimumAppliedCommand','MaximumAppliedCommand', ...
    'MinimumRawCommand','MaximumRawRate','MaximumAppliedRate','RateRMS','HJuseFraction', ...
    'FallbackDomain','FallbackResidual','PeakTorque'};
fid = fopen(file, 'w');
if fid < 0
    warning('HJ:Output', 'Cannot write metrics file: %s', file);
    return
end
fprintf(fid, '%s\n', strjoin(names, ','));
for i = 1:numel(cases)
    a = cases{i};
    if isempty(a), continue; end
    m = a.metrics;
    vals = [a.Bd, a.lambda_h, a.Rnu, m.thetaRMSE, m.tailThetaRMSE, m.finalStateErrorNorm, ...
        m.finalPhysicalStateErrorNorm, m.finalAppliedCommandError, m.finalRawCommandError, ...
        m.finalPressureError, m.finalHysteresisError, m.settlingTime, m.stateCost, m.rateCost, ...
        m.finiteHorizonCost, m.controlMin, m.controlMax, m.rawCommandMin, m.maximumRate, ...
        m.maximumAppliedRate, m.rateRMS, m.hjUseFraction, m.fallbackDomain, m.fallbackResidual, ...
        m.peakDisturbanceTorque];
    fprintf(fid, '"%s","%s","%s"', a.name, a.controller, a.status);
    fprintf(fid, ',%.16g', vals);
    fprintf(fid, '\n');
end
fclose(fid);
end
