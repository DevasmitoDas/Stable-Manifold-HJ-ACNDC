function a = f_smooth_abs(x, e)
% ==================================================================
% a = f_smooth_abs(x, e)
%
% F_SMOOTH_ABS  smooth surrogate of |x| used in the Bouc-Wen law (paper eq. (hyst))
%     s_e(x) = sqrt(x^2 + e^2) - e,   s_e(0) = 0,   0 <= |x| - s_e(x) <= e
%   evaluated in the rationalised form x^2/(sqrt(x^2+e^2)+e), which avoids the
%   cancellation near x = 0 (the same form as in symbSetting.m).
% ==================================================================
a = x.^2./(sqrt(x.^2 + e^2) + e);
end
