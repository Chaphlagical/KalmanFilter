function Minv = SPDinvIP(M);
%
% Grewal & Andrews, Kalman Filtering Theory and Practice Using MATLAB, 4th
% Edition, Wiley, 2014.
%
% Inverts symmetric positive definite matrix M IN PLACE using modified
% Cholesky decomposition M = U*D*U'
%
% INPUT: M, a symmetric positive definite matrix.
%
% OUTPUT: Minv, the inverse of M
%
% If the output argument (Minv) equals the input argument (M), inversion
% is done in place, overwriting M with its inverse.  The input matrix is
% destroyed in all cases.
%
[m,whatever] = size(M);
M = UD_decomp(M);  % modified Cholesky decomposition in-place
M = UDinv(M);     % inversion of modified Cholesky factors in-place
%
% Composing inv(M) from inv(D) and inv(U)
%
for j=m:-1:1,
    for i=j:-1:1,
        if i==j
            s = M(i,i);
        else
            s = M(i,i)*M(i,j);
        end;
        for k=1:i-1,
            s = s + M(k,i)*M(k,k)*M(k,j);
        end;
        Minv(i,j) = s;
        Minv(j,i) = s;
    end;
end;