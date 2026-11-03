function E = f_poly_exponents(degree, nx)
% ==================================================================
% E = f_poly_exponents(degree, nx)
%
% Exponents alpha of all monomials y^alpha with 1 <= |alpha| <= degree,
% one monomial per ROW (nx columns). For degree 3 and nx = 5 there are
% 5 + 15 + 35 = 55 monomials (paper eq. (ansatz)). The constant term is
% excluded, so g(0) = 0 for every coefficient vector.
% ==================================================================
if nargin < 2, nx = 5; end
E = zeros(0, nx);
for total = 1:degree
    E = [E; exponent_compositions(total, nx)]; %#ok<AGROW>
end
end

function E = exponent_compositions(total, nx)
if nx == 1
    E = total;
    return
end
E = zeros(0, nx);
for a = 0:total
    tail = exponent_compositions(total - a, nx - 1);
    E = [E; repmat(a, size(tail,1), 1) tail]; %#ok<AGROW>
end
end
