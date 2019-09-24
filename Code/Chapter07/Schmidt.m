function CPout = Schmidt(PhiCPin,GCQ)
% function CPout = Schmidt(PhiCPin,GCQ)
%
% Schmidt-Householder temporal update of triangular Cholesky factor of the
% covariance matrix of estimation uncertainty.
%
%   M. S. Grewal and A. P. Andrews,
%   Kalman Filtering: Theory and Practice Using MATLAB,
%   Wiley, 2014
%
% This algorithm performs Householder upper-triangularization of the matrix 
%
%               [ Phi*CPin, G*CQ ]
%
% where
%
%       Phi     is the n-by-n state transition matrix Phi[k];
%       CPin    is an n-by-n Cholesky factor of P[k-1+], the n-by-n
%               covariance matrix of a posteriori state estimation
%               uncertainty;
%       G       is the n-by-r process noise distribution matrix;
%       CQ      is an r-by-r Cholesky factor of Q[k-1], the r-by-r
%               covariance matrix of process noise.
%
% The inputs are the matrix products
%
%       PhiCPin = Phi*CPin
%       GCQ     = G*CQ
%
% The resulting output is
%
%       CPout,  the n-by-n upper triangular Cholesky factor of P[k-], the
%               covariance of a priori state estimation uncertainty.
%
% This method for the temporal update of Cholesky factors of covariance
% matrices was first implemented by Stanley F. Schmidt.  Prof. Thomas
% Kailath at Stanford University may have been the first to propose it.
%
[n,r] = size(GCQ);
for k = n:-1:1,
    sigma = 0; 
    for j = 1:r, 
        sigma = sigma + GCQ(k,j)^2; 
    end; 
    for j = 1:k, 
        sigma = sigma + PhiCPin(k,j)^2; 
    end; 
    alpha = sqrt(sigma); 
    sigma = 0; 
    for j = 1:r, 
        w(j)  = GCQ(k,j); 
        sigma = sigma + w(j)^2; 
    end; 
    for j = 1:k, 
        if j == k 
             v(j) = PhiCPin(k,j) - alpha; 
        else 
             v(j) = PhiCPin(k,j); 
        end; 
        sigma = sigma + v(j)^2; 
    end; 
    alpha = 2/sigma; 
    for i = 1:k, 
        sigma = 0; 
        for j = 1:r, 
            sigma = sigma + GCQ(i,j)*w(j); 
        end; 
        for j = 1:k, 
            sigma = sigma + PhiCPin(i,j)*v(j); 
        end; 
        beta = alpha*sigma; 
        for j = 1:r, 
            GCQ(i,j) = GCQ(i,j) - beta*w(j); 
        end; 
        for j = 1:k, 
             PhiCPin(i,j) = PhiCPin(i,j) - beta*v(j); 
        end; 
    end; 
end;
CPout = PhiCPin;
return;