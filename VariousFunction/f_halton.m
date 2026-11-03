function H = f_halton(n, bases, offset)
% ==================================================================
% H = f_halton(n, bases, offset)
%
% F_HALTON  Halton points (Halton 1960) in [0,1)^d, one point per column
%   n      : number of points
%   bases  : prime bases, one per coordinate (e.g. [2 3 5 7 11])
%   offset : the k-th point uses the index k + offset (leaped start)
% Used for the seed directions (Decide_xi.m), the collocation set and the
% held-out validation states of the value fit.
% ==================================================================
H = zeros(numel(bases), n);
for j = 1:numel(bases)
    for k = 1:n
        index = k + offset; fraction = 1/bases(j); value = 0;
        while index > 0
            value = value + fraction*mod(index, bases(j));
            index = floor(index/bases(j)); fraction = fraction/bases(j);
        end
        H(j,k) = value;
    end
end
end
