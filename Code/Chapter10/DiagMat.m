function D = DiagMat(v)
%
% Computes diagonal matrix with diagonal values given as a vector
%
% INPUT
%   v   a vector
% OUTPUT
%   D   a diagonal matrix
%
L   = length(v);
D   = zeros(L);
for k=1:L,
    D(k,k)  = v(k);
end;