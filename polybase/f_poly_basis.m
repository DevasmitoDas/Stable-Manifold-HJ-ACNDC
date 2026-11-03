function [phi, dphi] = f_poly_basis(Y, E)
% ==================================================================
% [phi, dphi] = f_poly_basis(Y, E)
%
% Monomials phi(b,i) = prod_j Y(j,i)^E(b,j) and their derivatives
% dphi(j,b,i) = d phi(b,i) / d Y(j,i).
%   Y : nx-by-np states (columns),  E : nb-by-nx exponents (f_poly_exponents)
% ==================================================================
nb = size(E,1); np = size(Y,2); nx = size(Y,1);
phi = ones(nb, np); dphi = zeros(nx, nb, np);
for j = 1:nx
    phi = phi.*(Y(j,:).^E(:,j));
end
for j = 1:nx
    ep = E; ep(:,j) = max(0, ep(:,j) - 1); v = ones(nb, np);
    for k = 1:nx
        v = v.*(Y(k,:).^ep(:,k));
    end
    v = v.*E(:,j);
    dphi(j,:,:) = reshape(v, 1, nb, np);
end
end
