function C = setdiff0(A, B, tol)

%SETDIFF0  Set difference with tolerance
% function x = solve0(varargin)
%
% (C) Michiel Hochstenbach 2021

if nargin < 3 || isempty(tol), tol = 1e-6; end
C = [];
for i = 1:length(A)
  in = true;
  for j = 1:length(B)
    if (A(i)==0 && B(j)==0) || abs(1-A(i)/B(j)) < tol, in = false; break; end; end
  if in, C(end+1) = A(i); end
end
