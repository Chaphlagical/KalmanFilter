function [yhat,Pyy] = UTsimplex(xhat,Cxx,trans,Wzero,alpha)
%
% 
% M. S. Grewal & A. P. Andrews
% Kalman Filtering Theory and Practice Using MATLAB
% 4th Edition, Wiley, 2014.
% 
% This function implements the "scaled simplex" unscented transform
% of Julier and Uhlmann.
%
% Modified after the Scilab code in Appendix A of S. J. Julier,
% "Reduced sigma point filters for the propagation of means and covariances
% through nonlinear transformations," available online at
% citeseer.ist.psu.edu/240235.html
%
% Calling sequence:
%
%   [yhat,Pyy] = UTsimplex(xhat,C,@trans,Wzero,alpha);
%
% Inputs:
%	xhat    Input mean vector
%   Cxx     Cholesky factor of covariance matrix
%   trans   Transformation function (passed as "@trans")
%	Wzero   Value of W(0), the weight applied to the mean (default = .5)
%	alpha   Scaling parameter (default = 1, the unscaled case)
%
% Outputs:
%	yhat    mean of y = trans(x)
%	Pyy     Covariance of y=trans(x).
%
%   The n+2 sample points will be
%
%       chi(0)      = hat(x), the mean
%       chi(1:n+1)  = hat(x) + columns of Chi, where
%       Chi         = Cxx*Xpts
%       Cxx*Cxx'    = Pxx, the state space covariance matrix
%
% Set default alpha and W(0), if necessary 
%
if (nargin < 3)
    error('Insufficient number of inputs');
elseif (nargin == 3)
    Wzero = .5;
    alpha = 1;
elseif (nargin == 4)
    alpha = 1;
end;
%
% Dimension the intermediate arrays
%
n           = length(xhat);
noPoints    = n+2;               % number of samples
Wpts        = zeros(1,noPoints); % sample weightings
Xpts        = zeros(n,noPoints); % samples
%
% Calculate the weight vector
%
Wpts(1)             = 1+(Wzero-1)/alpha^2;
Wpts([3:noPoints])  = (2.^[-n:-1]*(1-Wpts(1))); 
Wpts(2)             = Wpts(3);
%
% Calculate intermediate sigma points
%
for i=1:n,
    x           = 1/sqrt(2^i*Wpts(2)); 
    for j=2:i+1,
        Xpts(i,j)   = -x;
    end;
    Xpts(i,i+2) = x;
end;
%
% Calculate the actual sample points
%
Xpts = xhat*ones(1,n+2) + Cxx*Xpts;
%
% Transform the sample points
%
Ypts = trans(Xpts);
m    = size(Ypts,1);
%
% Calculate the output mean and covariance
%
yhat = zeros(m,1);
for k=1:n+2,
    yhat = yhat + Wpts(k)*Ypts(:,k);
end;
Pyy = zeros(m);
for k=1:n+2,
    Pyy = Pyy + Wpts(k)*(Ypts(:,k)-yhat)*(Ypts(:,k)-yhat)';
end;
return
