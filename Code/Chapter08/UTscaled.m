function [yhat,Pyy] = UTscaled(xhat,Cxx,trans,kappa,alpha,beta)
%
% 
% M. S. Grewal & A. P. Andrews
% Kalman Filtering Theory and Practice Using MATLAB
% 4th Edition, Wiley, 2014
% 
% This function implements the "scaled" unscented transform of Julier and
% Uhlmann, which approximates the mean and covariance of a variate
% y = trans(x), where "trans" is a nonlinear transformation and x is a
% variate with mean xhat and covariance Pxx = Cxx*Cxx'
%
% Calling sequence:
%
%   [yhat,Pyy] = UTscaled(xhat,C,@trans,kappa,alpha,beta);
%
% Inputs:
%	xhat    Input mean vector
%   Cxx     Cholesky factor of covariance matrix
%   trans   Transformation function (passed as "@trans")
%	kappa   Mean weight scaling (default = 0)
%	alpha   Sample spread parameter (default = 1, the unscaled case)
%	beta    Additional covariance scaling (default = 2 for Gaussian priors)
%
% Outputs:
%	yhat    mean of y = trans(x)
%	Pyy     Covariance of y=trans(x).
%
% Set default scaling parameters kappa, alpha and beta --- if necessary.
%
if (nargin < 3)
    error('Insufficient number of inputs');
elseif (nargin == 3)
    kappa = 0;
    alpha = 1;
    beta  = 2;
elseif (nargin == 4)
    alpha = 1;
    beta  = 2;
elseif (nargin == 5)
    beta  = 2;
end;
%
% Dimension the weight and sample arrays
%
n           = length(xhat);
noPoints    = 2*n+1;             % number of samples
Wpts        = zeros(1,noPoints); % sample weightings
Xpts        = zeros(n,noPoints); % samples
%
% Calculate the weight vector
%
    lambda          = alpha^2*(n + kappa) - n;
Wpts(1)             = lambda/(n+kappa);
for i=1:2*n,
    Wpts(i+1)       = 1/(2*(n+lambda));
end;
%
% Calculate sigma points
%
c                   = sqrt(n+lambda);
Xpts(:,1)           = xhat;
Xpts(:,2:n+1)       = xhat*ones(1,n)+c*Cxx; 
Xpts(:,n+2:2*n+1)   = xhat*ones(1,n)-c*Cxx; 
%
% Transform the sample points
%
Ypts = trans(Xpts);
m    = size(Ypts,1);
%
% Calculate the output mean and covariance
%
yhat    = (Wpts(1)+(1-alpha^2+beta))*Ypts(:,1)+Ypts(:,2:2*n+1)*Wpts(2:2*n+1)';
Pyy     = zeros(m);
for k=1:2*n+1,
    Pyy = Pyy + Wpts(k)*(Ypts(:,k)-yhat)*(Ypts(:,k)-yhat)';
end;
return
