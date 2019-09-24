function SensorOutput = hCompSens(x)
%
% function SensorOutput = hCompSens(x)
%
% Computes noise-free sensor output as a function of the state vector
%
%         [ scale factor ]
%     x = [ bias         ]
%         [ input        ]
%
%   Grewal & Andrews, Kalman Filtering: Theory and Practice Using MATLAB,
%   4th edition, Wiley, 2014.
%
SensorOutput = x(1)*x(3) + x(2);