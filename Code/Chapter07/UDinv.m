function InvUD = UDinv(UD)
%
% Grewal & Andrews, Kalman Filtering Theory and Practice Using MATLAB, 3rd
% Edition, Wiley, 2008.
%
% Inversion of UD decomposition factors M = U*D*U' of a symmetric
% positive definite matrix M, where
%
%   U is a unit upper triangular matrix
%     (i.e., ones on diagonal, zeros below)
%   D is a diagonal matrix with positive diagonal elements
%
%    inv(M) = inv(U*D*U') = U'\D\inv(U)
%
% INPUT: UD, a square matrix with the diagonal elements of D and
%        off-diagonal elements of U.
%
% OUTPUT: InvUD, a square matrix with the diagonal elements of inverse(D)
%         and off-diagonal elements of Inverse(U)
%
% If the output argument (InvUD) equals the input argument (UD), inversion
% is performed in-place, overwriting UD with InvUD.
%
[m,whatever] = size(UD);
for i=m:-1:1,
    InvUD(i,i) = 1/UD(i,i);
    for j=m:-1:i+1, 
        InvUD(i,j) = -UD(i,j); 
        for k=i+1:j-1, 
            InvUD(i,j) = InvUD(i,j) - UD(i,k)*InvUD(k,j); 
        end;
        UDinv(j,i) = 0;
    end; 
end;
