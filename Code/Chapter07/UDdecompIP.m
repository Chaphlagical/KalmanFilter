function M = UDdecompIP(M)
%
% Performs UD decomposition of a symmetric positive definite matrix in
% place.
%
% INPUT: Symmetric positive definite m-by-m matrix M.
% OUTPUT: M is overwritten with U (above the diagonal) and D (below the
% diagonal) such that M = U*D*U'
%
[rows,cols] = size(M);
if (rows ~= cols)|(rows < 1)|(nargin ~= 1)
    error('Input to UDdecompIP must be a square matrix.');
else
    m = rows;
    for j = m : -1:1, 
        for i = j : -1:1, 
            sigma = M (i, j); 
            for k = j + 1 : m, 
                sigma = sigma - M(i, k) *M(k, k) *M(j, k); 
            end; 
            if i == j 
                M(j, j) = sigma; 
            else 
                M(i, j) = sigma/M(j, j); 
            end; 
        end; 
    end;
end;