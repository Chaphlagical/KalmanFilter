function [xsk,Psk] = BMFLS(Lag,k,zk,xskminus1,Pskminus1,Phi,Q,H,R)
%
% Grewal & Andrews,
% Kalman Filtering: Theory and Practice Using MATLAB, 4th Edition,
% Wiley, 2014
%
% Biswas-Mahalanabis Fixed-Lag Smoother with inputs
% 
%   Lag       = maximum number of lag time-steps (must be >= 0)
%   k         = measurement number from start (=1 at start),
%   zk        = measurement vector (column),
%   xskminus1 = previous a posteriori smoother estimate (column vector),
%   Pskminus1 = covariance of same,
%   Phi       = state transition matrix of filter,
%   Q         = dynamic noise covariance of filter,
%   H         = measurement sensitivity matrix of filter,
%   R         = measurement noise covariance of filter
%
% and outputs
%
%   xsk       = current a posteriori smoother estimate (column vector),
%   Psk       = covariance of same,
%
% Determine
%           m = number of measurement variables.
%           n = number of  state variables
%
[m,n] = size(H);
%
% Compute the cutoff block index for state vector augmentation
%
if k <= Lag+1
    cutoff = k-1;
else
    cutoff = Lag;
end;
%
% Temporal update
%
xsk   = [Phi*xskminus1(1:n);xskminus1(1:cutoff*n)];
Psk   = [Phi*Pskminus1(1:n,1:n)*Phi'+Q,Phi*Pskminus1(1:n,1:n*cutoff);Pskminus1(1:n,1:n*cutoff)'*Phi',Pskminus1(1:n*cutoff,1:n*cutoff)];
%
%  Compute the Kalman gain
%
PHT    = Psk(:,1:n)*H';
Rprime = H*PHT(1:n,1:m)+R;
K      = PHT/Rprime;
%
% Expected measurement and innovations
%
zhat  = H*xsk(1:n);
nu    = zk - zhat; % innovations
%
% Observational update
%
xsk  = xsk + K*nu;
Psk  = Psk - K*PHT';
Psk  = 0.5*(Psk+Psk');
return;

