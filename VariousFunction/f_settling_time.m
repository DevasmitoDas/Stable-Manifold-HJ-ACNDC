function ts = f_settling_time(t, e, band)
% ==================================================================
% ts = f_settling_time(t, e, band)
%
% F_SETTLING_TIME  instant after the last exit from |e| <= band
%   (0 if e never leaves the band, Inf if it is outside at the final time)
% ==================================================================
last = find(abs(e) > band, 1, 'last');
if isempty(last)
    ts = 0;
elseif last == numel(t)
    ts = inf;
else
    ts = t(last+1);
end
end
