function [x1,P11,P12,P22] = SKF_Temp(x1,P11,P12,P22,Phi1,Phi2,Q1,Q2);
% 
% Schmidt-Kalman Filter Temporal Update.
% 
% M. S. Grewal & A. P. Andrews
% Kalman Filtering: Theory and Practice Using MATLAB, 4th Edition,
% Wiley, 2014.
%
% INPUTS
%    x1     a posterior estimate of state subvector to be estimated
%    P11    a posteriori covariance of estimated state sub vector uncertainty
%    P22    a posteriori covariance of non-estimated state vector uncertainty x2
%    P12    cross-covariance between estimated and non-estimated components
%    Phi1  state transition submatrix for estimated state subvector x1
%    Phi2  state transition submatrix for non-estimated state subvector x2
%    Q1    dynamic disturbance covariance for estimated subvector
%    Q2    dynamic disturbance covariance for non-estimated subvector
%
% OUTPUTS
%    x1     a priori estimate of state subvector to be estimated
%    P11    a priori covariance of estimated state sub vector uncertainty
%    P22    a priori covariance of non-estimated state vector uncertainty x2
%    P12    a priori cross-covariance 
%
% Schmidt-Kalman temporal update implementation
%
x1    = Phi1*x1;
P11   = Phi1*P11*Phi1' + Q1;
P12   = Phi1*P12*Phi2';
P22   = Phi2*P22*Phi2' + Q2;
%
% Maintain symmetry of diagonal blocks of covariance matrix
%
P11   = .5*(P11+P11');        
P22   = .5*(P22+P22');
