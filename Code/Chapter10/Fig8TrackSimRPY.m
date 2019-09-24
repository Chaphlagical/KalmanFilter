function [xENU,vENU,aSensedRPY,omegaSensedRPY,CENU2RPY,Roll,Pitch,Yaw] = Fig8TrackSimRPY(t)
%
% Simulated dynamic conditions on a figure-8 track of 
%
% Track length = 10^4/6 meters = 10/6 km
% Crossover bridge height      = 20 meters
% Lap time                     = 60 seconds
%
% INPUT
%
%   t = time from undercrossing
%
% OUTPUTS (all column vectors representing INS state in SI units)
%
%   xENU        = [Easting;Northing;Altitude] with respect to undercrossing location
%   vENU        = [vE;vN;vU] ENU velocity with respect to track
%   aSensedRPY  = sensed acceleration in RPY coordinates [m/s/s]
%   omegaSensedRPY  = sensed rotation rates in RPY coordinates [rad/s]
%   CENU2RPY    = coordinate tranform matrix from ENU to RPY coordinates
%
%
% PARAMETER SETTINGS
%
S            = 107.63556307;    % [m] scaling for 10/6 km track length
Tlap         = 60;              % [s]
BridgeHeight = 20;              % [m]
MaxBank      = 0.467765692827;  % [rad] ~ 26.801 [deg]
LatDeg       = 39.7931;         % Approximate location of 
LonDeg       = -86.2389;        % Indianapolis Motor Speedway
LatRad       = LatDeg*pi/180;
LonRad       = LonDeg*pi/180;
sLat         = sin(LatRad);
cLat         = cos(LatRad);
omegaEarth   = 0.7292115e-4 * [0;cLat;sLat]; % Earthrate in ENU coordinates
%
s = sin(0.2e1*pi*t/Tlap);
c = cos(0.2e1*pi*t/Tlap);
%
% Cartesian variables
%
xENU        = [0.2e1*S*s*c;0.3e1*S*s;BridgeHeight*(0.1e1 - c)/0.2e1];
vENU        = [0.4e1*S*pi*(c - s)*(c + s)/Tlap;0.6e1*S*c*pi/Tlap;BridgeHeight*s*pi/Tlap];
aENU        = [-0.32e2*S*pi^2*c*s/Tlap^2;-0.12e2*S*s*pi^2/Tlap^2;0.2e1*BridgeHeight*c*pi^2/Tlap^2];
speed       = norm(vENU);
speedsq     = speed^2;
omegaENU    = cross(vENU,aENU)/speedsq;
aSensedENU  = aENU + [0;0;9.8];
omegaSensedENU = omegaENU + omegaEarth;
%
%  Euler angles
%
heading = atan2(vENU(1),vENU(2));
Yaw     = pi/2 - heading;
Pitch   = asin(vENU(3)/speed);
Roll    = -MaxBank*s;
%
% RPY axis directions in ENU coordinates
%
CENU2RPY    = ENU2RPY(Roll,Pitch,Yaw);
aRPY        = CENU2RPY*aENU;              % RPY acceleraton wrt track
aSensedRPY  = CENU2RPY*aSensedENU;        % RPY sensed acceleration
omegaRPY    = CENU2RPY*omegaENU;          % RPY rotation rate wrt track
omegaSensedRPY = CENU2RPY*omegaSensedENU; % RPY sensed rotation rate
