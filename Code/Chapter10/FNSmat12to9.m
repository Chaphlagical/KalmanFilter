function FNS = FNSmat12to9(Csens2nav,a,omega)
% 
% M. S, Grewal & A. P. Andrews
% Kalman Filtering: Theory and Practice Using MATLAB, 4th edition
% Wiley, 2014
%
% Dynamic coupling matrix from 12-state sensor compensation parameter
% errors into 9-state inertial navigation errors.
% 
% The 12 sensor compensation parameters in question  are
%   3 accelerometer biases
%   3 accelerometer scale factors
%   3 gyro biases
%   3 gyro scale factors
% The 9 navigation error components are
%   3 position error components in navigation coordinates
%   3 velocity error components in navigation coordinates
%   3 rotation error components in navigation coordinates
%
% INPUTS:
%   Csens2nav   coordinate transformation matrix from sensor coordinates to
%               navigation coordinates (can be used for strapdown INS)
%   a           sensed accelerations in sensor coordinates
%   omega       sensed rotation rates in sensor coordinates
% 
% OUTPUTS:
%
%   FNS         9x12 dynamic coefficient sub-matrix coupling sensor
%               compensation errors into navigation errors
%
FNS = [zeros(3,12);
    Csens2nav,Csens2nav*DiagMat(a),zeros(3,6);
    zeros(3,6),Csens2nav,Csens2nav*DiagMat(omega)];
return;