function [K,Pout] = josephdv(z,R,H,P)
%
% Joseph "stabilized" Kalman filter measurement
% update as implemented by Thomas De Vries.
% 
% M. S. Grewal & A. P. Andrews
% Kalman Filtering Theory and Practice Using MATLAB
% Third Edition, Wiley & Sons, 2008
% 
%
T1 = P*H';
T2 = H*T1 + R;
K  = T1/T2;
T3 = .5*K*T2 - T1;
T4 = T3*K';
Pout = P + T4 + T4';

