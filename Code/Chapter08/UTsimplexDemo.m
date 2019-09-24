%
% Simplex unscented transform demonstration for the nonlinear function
%
%    f(x) = [x(4)*x(2);x(1)*x(3);x(2)*x(4)] % (function NonLin)
%
% using four alternative sample-sets from the initial covariance matrix P:
%
%  1. Eigenvectors of P scaled by square roots of eigenvalues
%  2. Columns of upper triangular Cholesky factor of P
%  3. Columns of lower triangular Cholesky factor of P
%  4. Columns of symmetric Cholesky factor of P
%
%  This should demonstrate superiority of UT over EKF, and the relatively
%  high sensitivity of the UT-estimated mean and covariance to the choice
%  of Cholesky factors used for the UT sampling.
%
% 
% M. S. Grewal & A. P. Andrews
% Kalman Filtering Theory and Practice Using MATLAB
% 4th Edition, Wiley, 2014.
% 
close all;
clear all;
%
disp('Demonstration of simplex unscented transform of nonlinear measurement');
disp('using different Cholesky factors for unscented transform.');
disp(' ');
%
disp('4x4 a priori covariance matrix of the state variable:');
%
T               = 1/2*[1,1,1,1;1,1,-1,-1;1,-1,1,-1;1,-1,-1,1]; % orthogonal
P               = T*[4,0,0,0;0,36,0,0;0,0,100,0;0,0,0,196]*T',
%
% Eigenvector Choleksy factor
%
[U,Lambda,V]    = svd(P);
CEV             = U*sqrt(Lambda);
%
% Symmetric Choleksy factor
%
Csym            = U*sqrt(Lambda)*U';
%
% Upper triangular Choleksy factor
%
CUT             = utchol(P);
%
% Lower triangular Choleksy factor
%
CLT             = chol(P)';
%
% Initial mean
%
disp('A priori mean of state variable:');
xhat = [0;0;0;0],
%
% True measurement mean from NonLin transformation of xhat(1:3), P(1:3,1:3)
%
Tzhat    = zeros(3,1);
Tzhat(1) = P(2,3)*(1 - P(2,3)^2 / P(2,2) / P(3,3))^(-3/2);
Tzhat(2) = P(3,1)*(1 - P(3,1)^2 / P(3,3) / P(1,1))^(-3/2);
Tzhat(3) = P(1,2)*(1 - P(1,2)^2 / P(1,1) / P(2,2))^(-3/2);
disp('True mean of nonlinear measurement:')
Tzhat,
%
% Extended Kalman filter transformation of mean and variance
%
disp('Measurement mean and covariance using extended Kalman filter:');
Lzhat = NonLin(xhat(1:3)),
H     = [0,xhat(3),xhat(2),0;xhat(3),0,xhat(1),0;xhat(2),xhat(1),0,0];
LPzz  = H*P*H',
%
% simplex unscented transform solutions
%
disp('UT mean and covariance using scaled eigenvector samples:')
[sUTzhatEV,sUTPzzEV]    = UTsimplex(xhat(1:3),CEV(1:3,1:3),@NonLin,.5,1),
disp('UT mean and covariance using upper triangular Cholesky factors:')
[sUTzhatUT,sUTPzzUT]    = UTsimplex(xhat(1:3),CUT(1:3,1:3),@NonLin,.5,1),
disp('UT mean and covariance using lower triangular Cholesky factors:')
[sUTzhatLT,sUTPzzLT]    = UTsimplex(xhat(1:3),CLT(1:3,1:3),@NonLin,.5,1),
disp('UT mean and covariance using symmetric Cholesky factors:')
[sUTzhatsy,sUTPzzsy]    = UTsimplex(xhat(1:3),Csym(1:3,1:3),@NonLin,.5,1),
