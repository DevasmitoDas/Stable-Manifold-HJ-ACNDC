% ==================================================================
% P3: x(t), p(t) を求めるための初期値ξの決定 / seeds xi
%
% The seeds live in the DIAGONALISED coordinates s (stable block of
% (s,chi) = T^{-1}(y,mu)); a seed is therefore not a prescribed initial state:
% y(0) = xi + L*chi(0) (paper Sec. 3.2).
%   iniFlag = 1 : Decide_xi.m in the problem folder (+-e_i and 8 Halton
%                 directions of radius radi -> 18 seeds), as in the paper.
% The reference sphere sampling (iniFlag = 0, theta1, theta2, ...) is not used
% for this problem.
% ==================================================================

if iniFlag == 1
	Decide_xi       % defines xi_dir (unit directions) and xi_sum (seeds, one per row)
else
	error('HJ:P3', 'Set iniFlag = 1 in systemData.m: this problem uses Decide_xi.m.');
end

[ini_row, ini_col] = size(xi_sum);
total = ini_row*(knum + 1);
fprintf('\nNumber of seeds xi : %d (radius %g)\n', ini_row, radi)

% -------------------------- report ------------- %
Report_P3 = {'num of xi_sum', ini_row, 'paper: 18'; ...
	'radius of xi', radi, 'paper: 2e-4'; ...
	'det of Trs (eig vectors of F)', det(Trs), 'be careful if 0'};
