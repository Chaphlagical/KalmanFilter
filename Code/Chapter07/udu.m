function [U,D] = udu(P)
%
% Grewal & Andrews, Kalman Filtering Theory and Practice Using MATLAB, 4th
% Edition, Wiley, 2014.
%
% function [U,D] = udu(P)
%
% Performs modified Cholesky decomposition of symmetric positive-definite
% matrix P (input).
%
% Returns the modified Cholesky factors U and D such that
%
%   1. U*D*U' == P
%   2. U is unit upper triangular, i.e., U is a square matrix of the same
%      dimensions as P, with zeros below its main diagonal and ones on its
%      main diagonal.
%   3. D is a diagonal matrix, i.e., a square matrix of the same dimensions
%      as P, with zeros above and below its main diagonal.
%
[rows,cols] = size(P);
if rows == cols
    m = rows;
    P = (P+P')/2;
    for j = m:-1:1, 
        for i = j:-1:1, 
            sigma = P(i,j); 
            for k = j+1:m, 
                sigma = sigma - U(i,k)*D(k,k)*U(j,k); 
            end; 
            if i == j 
                D(j,j) = sigma; 
                U(j,j) = 1; 
            else 
            U(i,j) = sigma/D(j,j); 
            end; 
        end; 
    end; 
else
    error('Input matrix to function udu must be square');
end;