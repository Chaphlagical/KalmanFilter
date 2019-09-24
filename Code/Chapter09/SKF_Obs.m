function [x1,P11,P12,P22] = SKF_Obs(x1,P11,P12,P22,z,H1,H2,R);
% 
% Schmidt-Kalman Filter Observational Update.
% 
% M. S. Grewal & A. P. Andrews
% Kalman Filtering: Theory and Practice Using MATLAB, 4th Edition, 
% Wiley, 2014.
%
% INPUTS
%    x1   a priori estimate of state subvector to be estimated
%    P11  a priori covariance of estimated state subvector uncertainty x1
%    P22  a priori covariance of non-estimated state subvector uncertainty x2
%    P12  a priori cross-covariance between x1 and x2 estimation errors
%    z    measurement
%    H1   measurement sensitivity to estimated state subvector x1
%    H2   measurement sensitivity to non-estimated state subvector x2
%    R    covariance of measurement noise
%    
% OUTPUTS
%    x1   a posteriori estimate of state subvector to be estimated
%    P11  a posteriori covariance of estimated state subvector uncertainty x1
%    P22  a posteriori covariance of non-estimated state subvector uncertainty x2
%    P12  a posteriori cross-covariance between x1 and x2 estimation errors
%
% Calculate Schmidt-Kalman gain
%
C0  = P11*H1'+P12*H2';
D0  = P12'*H1'+P22*H2';
E0  = H1*C0+H2*D0+R;
Ksk = C0/E0;
%
% Schmidt-Kalman observational update
%
x1    = x1 + Ksk*(z - H1*x1);
A0    = eye(2) - Ksk*H1;
B0    = Ksk*H2;
P12o  = A0*P12 - B0*P22;
P11o  = (A0*P11 - B0*P12')*A0' - P12*B0' + Ksk*R*Ksk';
P11   = .5*(P11o+P11o');
P12   = P12o;