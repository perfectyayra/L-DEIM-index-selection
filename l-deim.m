function [C,R,M] = ldeim(A,U,V,r, k)

%CUR_DEIM_SOFT  Soft DEIM incurred CUR decomposition
% function [C, R, M] = ldeim(A,U,V, r,k)
% A is the data matrix
% U contains the left singular vectors
% V contains the right singular vectors
% r is the desired rank of the approximation
% k the number of available singular vectors
% k << r
% A ~ CMR,  C = A(:,icol);  R = A(irow,:)
%
% (C) Perfect Gidisu, Michiel Hochstenbach 2021


irow=zeros(1,k);
icol=zeros(1,k);
for j = 1:k
  % Soft indexing
  [~, irow(j)] = max(abs(U(:,j)));
  [~, icol(j)] = max(abs(V(:,j)));
  if j<k
    U(:,j+1) = U(:,j+1) - U(:,1:j) * (U(irow(1:j),1:j) \ U(irow(1:j),j+1));
    V(:,j+1) = V(:,j+1) - V(:,1:j) * (V(icol(1:j),1:j) \ V(icol(1:j),j+1));
  end
end
U2 = rownorms(U); % Re-evaluate index set
V2= rownorms(V);
[~, irow2] = sort(U2, 'descend');
[~, icol2] = sort(V2, 'descend');
irow2 = setdiff0(irow2, irow); 
icol2 = setdiff0(icol2, icol); 
irow = [irow(1:k) irow2(1:r-k)];
icol = [icol(1:k) icol2(1:r-k)];
C= A(:,icol);
R= A(irow,:);
M = C \ (A / R);
