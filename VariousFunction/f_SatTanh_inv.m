function u = f_SatTanh_inv(ua, umin, umax)
% ==================================================================
% u = f_SatTanh_inv(ua, umin, umax)
%
% F_SATTANH_INV  inverse of the smooth saturation, defined for umin < ua < umax
%     u = u_m + u_h*atanh((ua - u_m)/u_h)
%   e.g. u* = sigma^{-1}(0.6) = 0.5 + 0.5*atanh(0.2) = 0.6013662770
% ==================================================================
um = (umin + umax)/2;  uh = (umax - umin)/2;
u = um + uh*atanh((ua - um)/uh);
end
