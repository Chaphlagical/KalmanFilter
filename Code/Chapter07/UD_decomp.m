function UD = UD_decomp(M)
%
% Grewal & Andrews, Kalman Filtering Theory and Practice Using MATLAB, 3rd
% Edition, Wiley, 2008.
%
% UD decoposition in-place
%
% INPUT: M, a symmetric positive definite real matrix
%
% OUTPUT: UD, an upper triangular matrix, the diagonal of which is the
%         diagonal of a square matrix D, and the strictly upper diagonal
%         of which is the strictly upper diagonal part of a unit upper
%         triangular matrix U such that M = U*D*U'
%
% If the output argument (UD) is M, decomposition is performed in-place.
%
[m,whatever] = size(M);
for j=m:-1:1,
    for i=j:-1:1,
        sigma = M(i,j);
        for k=j+1:m,
            sigma = sigma - UD(i,k)*UD(k,k)*UD(j,k);
        end;
        if i==j
            UD(j,j) = sigma;
        else
            UD(i,j) = sigma/UD(j,j);
            UD(j,i) = 0;
        end;
    end;
end;
            
        