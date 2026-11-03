function [Ak,T] = ordblk(A,cut,ordtype)

% Real ordered block diagonal form 
%
% [Ak,T] = ordblk(A,cut,ordtype) produces a real block-ordered decomposition
% of matrix A such that                             
%
%             Ak = blkdiag([Jordan Block of repeated zero eigenvalues],
%                          [complex eigenvalue]   , .....]);
%                                              2x2
%             where Jordan Block of the repeated zero eigenvalues is in the form of
%               
%                          [0 1 0 0 ...
%                           0 0 1 0 ...
%                           0 0 0 1 ...
%                            ..........]
%
%             and each complex eigenvalue a+bj is in 2x2 real matrix
%
%                          [a b;
%                          -b a].
%
% The resulting matrix Ak is ordered in block ascending form based on their 
% eigenvalue magnitudes.
% ORDTYPE = 'mag'  ( abs(eig(G1))  < abs(eig(G2)), default)
% ORDTYPE = 'real' ( real(eig(G1)) > real(eig(G2)) )

%
% See also SLOWFAST, STABPROJ, and BLKRSCH.

% R. Y. Chiang & M. G. Safonov 1/2003
% Copyright 1988-2004 The MathWorks, Inc. 
%       $Revision: 1.1.6.7 $
% All Rights Reserved.
nag1 = nargin;

if nag1 == 2
    ordtype = 'mag';
end

n = length(A);
cut = n; 

if nag1 < 2
   cut = n;
   ordtype = 'mag';
end 

od = ordtype(1:3);
switch od
 case 'mag'
    ordtype = 6; % ascending order
 case 'rea'
    ordtype = 3; % descending order
end

if cut == n
   if ordtype == 6
        [T,Ak] = reig(A,2);      % magnitude ascending order
   end
   if ordtype == 3
        [T,Ak] = reig(A,1);      
        T = T(:,n:-1:1);
        Ak = Ak(n:-1:1,n:-1:1);  % real part desecend order
   end
   return
end

%
% General case:
%
return
[Tk,Aks] = blkrsch(A,ordtype,cut);   

Ak1 = Aks(1:cut,1:cut); Ak2 = Aks(cut+1:n,cut+1:n); A12 = Aks(1:cut,cut+1:n);

X = lyap(Ak1,-Ak2,A12);
[Tk2,Ak2] = reig(Ak2,2); 
Tdiag = [eye(cut,cut) X*Tk2; zeros(n-cut,cut) Tk2];

Ak = blkdiag(Ak1,Ak2);
T = Tk*Tdiag;
%
% ---------- End of ORDBLK.M % RYC/MGS
