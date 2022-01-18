function n = rownorms(A, p)

%ROWNORMS  Row norms
% function n = rownorms(A, p)
%
% See also COLNORMS
%
% Revision date: February 7, 2021
% (C) Michiel Hochstenbach 2021

if nargin < 2 || isempty(p), p = 2; end

if p == 1
  n = sum(abs(A),2);
elseif p == 2
  n = zeros(size(A,1), 1);
  for j = 1:size(A,1), n(j) = norm(A(j,:)); end
else  % inf-norm
  n = max(abs(A),[],2);
end
