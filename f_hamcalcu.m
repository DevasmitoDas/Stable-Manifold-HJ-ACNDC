function [H, stage, term, nu, fs] = f_hamcalcu(Y, Mu, des)
% ==================================================================
% [H, stage, term, nu, fs] = f_hamcalcu(Y, Mu, des)
%
% F_HAMCALCU  minimised Hamiltonian (paper eq. (Hbar)) on column batches
%     H     = Hbar(y,mu) = mu'*a_sigma(y) + y'*Qy*y - (G_y'*mu)^2/(4*Rnu)
%     nu    = nu*(mu) = -G_y'*mu/(2*Rnu)               (eq. (nustar))
%     fs    = a_sigma(y) + G_y*nu
%     stage = y'*Qy*y + Rnu*nu^2   (running cost l),   term = mu'*fs
% H is evaluated from the substituted Hamiltonian of buildProblem.m.
% NOTE) columns are samples (the reference version used rows).
% ==================================================================
nu = des.prob.optimalRate(Y, Mu);
fs = des.prob.f(Y) + des.B*nu;
stage = sum(Y.*(des.Q*Y),1) + des.R*nu.^2;
term = sum(Mu.*fs,1);
H = des.prob.hamiltonian(Y, Mu);
end
