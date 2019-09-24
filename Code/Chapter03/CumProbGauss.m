%   M. S. Grewal and A. P. Andrews
%   Kalman Filtering Theory and Practice Using MATLAB, 4th Edition,
%   Wiley, 2014
%   Chapter 3
function cumP = CumProbGauss(x);
%
% Computes the cumulative probability function for the zero-mean
% unit-normal Gaussian distribution.
%
CumP = 0.5*(1 + erf(x/sqrt(2)));