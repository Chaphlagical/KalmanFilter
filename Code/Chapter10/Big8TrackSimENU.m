function [xENU,vENU,aENU,aSensedENU,omegaENU,omegaSensedENU] = Big8TrackSimENU(t,LatRad,LonRad,Alt)
% 
% M. S, Grewal & A. P. Andrews
% Kalman Filtering: Theory and Practice Using MATLAB, 4th edition
% Wiley, 2014
%
%
% Simulated East-North-Up dynamic conditions on a 100-km figure-8 track 
%
% INPUT
%
%   t       = time from undercrossing [s]
%   LatRad  = Latitude of track at crossover point
%   LonRad  = Longitude of track at crossover point
%   Alt     =  Altitude of track at bottom of crossover
%
%
% OUTPUTS (all column vectors representing INS state in SI units)
%
%   xENU            = [Easting;Northing;Altitude] with respect to undercrossing location
%   vENU            = [vE;vN;vU] velocity with respect to track
%   aSensedENU      = [aE;aN;aU] sensed acceleration with respect to track
%   omegaSensedENU  = sensed rotation rates, including earthrate
%                       [rad/s]
%
% PARAMETER SETTINGS
%
S            = 6458.1337842;      % [m] scaling for 100 km track length
Tlap         = 3600;              % [s] lap time (one hour)
BridgeHeight = 10;                % [m] crossover separation
MaxBank      = 0.00779609488045;  % [rad] ~0.44668 [deg]
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