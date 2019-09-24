function [Phi,Q,P] = MODL5Params(SigmaAcc,SigmaVel,TauAcc,DeltaT);
%
% Computes parameters of the MODL5 stochastic model, given desired values:
%
% INPUTS
%   SigmaAcc  = RMS acceleration
%   SigmaVel  = RMS velocity
%   TauAcc    = correlation time-constant of accelerations
%   DeltaT    = discrete time interval
%
% OUTPUTS   = Parameters of discrete-time model
%   Phi     = State transition matrix
%   Q       = Process noise covariance
%   P       = Steady-state covaraince matrix of state uncertainty
%
TauVel      = SigmaVel*(SigmaVel + sqrt(SigmaVel^2 + 4*TauAcc^2*SigmaAcc^2))/SigmaAcc^2/TauAcc/2;
RhoVelAcc   = SigmaVel/SigmaAcc/TauVel;
if abs(RhoVelAcc) >=1, error('MODL5Params: Incompatible inputs, |RhoVelAcc| >= 1'); end;
F           = [-1/TauVel, 1; 0, -1/TauAcc];
Phi         = expm(DeltaT*F);
P           = [SigmaVel^2, RhoVelAcc*SigmaVel*SigmaAcc; RhoVelAcc*SigmaVel*SigmaAcc,SigmaAcc^2];
Q           = P - Phi*P*Phi';
% 
% M. S, Grewal and A. P. Andrews,
% Kalman Filtering: Theory and Practice Using MATLAB, 4th Edition,
% Wiley, 2014
%

