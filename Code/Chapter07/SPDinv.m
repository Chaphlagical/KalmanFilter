function Minv = SPDinv(M);
%
% Grewal & Andrews, Kalman Filtering Theory and Practice Using MATLAB, 4th
% Edition, Wiley, 2014.
%
% Inverts symmetric positive definite matrix M using modified Cholesky
% decomposition M = U*D*U'
%
% INPUT: M, a symmetric positive definite matrix.
%
% OUTPUT: Minv, the inverse of M
%
% If the output argument (Minv) equals the input argument (M), inversion
% is done in place, overwriting M with its inverse.
%
[m,whatever] = size(M);
UD = UD_decomp(M);  % modified Cholesky decomposition
UD = UDinv(UD);     % inversion of modified Cholesky factors
%
% Composing inv(M) from inv(D) and inv(U)
%
for j=m:-1:1,
    for i=j:-1:1,
        if i==j
            s = UD(i,i);
        else
            s = UD(i,i)*UD(i,j);
        end;
        for k=1:i-1,
            s = s + UD(k,i)*UD(k,k)*UD(k,j);
        end;
        Minv(i,j) = s;
        Minv(j,i) = s;
    end;
end;