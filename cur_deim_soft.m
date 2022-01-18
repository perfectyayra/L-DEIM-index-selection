function irow = cur_deim_soft(U,r,k)

%CUR_DEIM_SOFT  Soft DEIM incurred CUR decomposition
% function [icol, irow, M] = cur_deim_soft(A, k, U, V, S)
%
% A ~ CMR,  C = A(:,icol);  R = A(irow,:)
%
% Revision date: February 6, 2021
% (C) Perfect Gidisu, Michiel Hochstenbach 2021



for j = 1:k
  % Soft indexing
  [~, irow(j)] = max(abs(U(:,j)));
  U(:,j+1:k) = U(:,j+1:k) - U(:,1:j) * (U(irow,1:j) \ U(irow,j+1:k));
end
U2 = rownorms(U(:,1:k)); % Re-evaluate index set
[~, irow2] = sort(U2, 'descend');
irow2 = setdiff0(irow2, irow); 
irow = [irow(1:k) irow2(1:r-k)];
%M = A(:,icol) \ (A / A(irow,:));
