function [Phik,Qk] = VanLoan(F,DeltaT,Qc,G)
%
%   M. S. Grewal and A. P. Andrews
%   Kalman Filtering Theory and Practice Using MATLAB, 4th Edition,
%   Wiley, 2014
%   Chapter 4
%
% function [Phki,Qk] = VanLoan(F,DeltaT,Qc,G)
%
% Given the time-invariant matrices F, Qc and G [optional] of the
% stochastic differential equation
%
% d/dt x(t) = F x + G w(t),
%
% the dynamic disturbance noise coraviance Qc = E<ww'>, and the time-step
% DeltaT for an equivalent process model in discrete time, finds the
% equivalent matrices Phik and Qk for the discrete-time implementation
%
% x[k] = Phik x[k-1] + w[k-1],
%
% where Qk is the covariance of the zero-mean white noise process {w[.]}.
%
[m,n] = size(F);
if m ~= n
    error('VanLoan input F must be a square matrix');
end;
if nargin == 4
    M    = DeltaT*[-F,G*Qc*G';zeros(n),F'];
elseif nargin == 3
    M    = DeltaT*[-F,Qc;zeros(n),F'];
else
    error('function VanLoan requires 3 or 4 arguments')
end
N    = expm(M);
Phik = N(n+1:2*n,n+1:2*n)';
Qk   = Phik*N(1:n,n+1:2*n);
return;