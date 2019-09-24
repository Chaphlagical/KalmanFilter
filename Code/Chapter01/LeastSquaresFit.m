close all;
clear all;
%
disp('GREWAL AND ANDREWS,');
disp('KALMAN FILTERING: THEORY AND PRACTICE USING MATLAB, 4TH EDITION,');
disp('WILEY, 2014.');
disp('CHAPTER 1');
disp('-------------------------------------------------------------------');
disp('LEAST-SQUARES FIT OF INERTIAL NAVIGATION RMS HORIZONTAL UNCERTAINTY');
disp('TO CEP RATE, WHERE CEP = "CIRCULAR ERROR PROBABLE" AND UNITS ARE');
disp('NAUTICAL MILES (1852 meters).');
disp('-------------------------------------------------------------------');
%
% 1. LOAD DATA (RHS horizontal error [m] from medium accuracy inertial
%               navigation system, solution from Riccati equation.)
%
load 'RMSHorINS1m.mat'
%
% CREATE TIMELINE
%
N       = length(RMSHorINS);
hr      = 4*(0:N-1)/(N-1);
plot(hr,RMSHorINS/1852,'k-','LineWidth',2); % CEP is "circle of equal probability,"
%                                             a circle around the estimate
%                                             of sufficient radius such
%                                             that it is equally likely
%                                             that the true value is inside
%                                             or outside of that circle.
%                                             Division by 1852 converts
%                                             from meters to nautical miles
hold on;
xlabel('TIME IN HOURS');
ylabel('RMS HORIZONTAL NAVIGATION ERROR IN NMi');
title('INPUT DATA  FROM CHAPTER 10 FOR LEAST-SQUARES STRAIGHT-LINE FIT');
%
% CALCULATE GRAMIAN (SCALAR)
%
G       = hr*hr';   % hr is a row vector.
%
% CALCULATE LEAST-SQUARES STRAIGHT-LINE SLOPE IN NMi/Hr
%
if G ~= 0
    disp(['GRAMIAN             = ',num2str(G),' Hr^2.']);
    CEPrate = hr*RMSHorINS'/G/1852;
    disp(['CALCULATED CEP RATE = ',num2str(CEPrate),' NMi/Hr.']);
    disp('(Shown by dashed line on the plot.)');
    disp('(THIS WILL MAKE MORE SENSE IN CHAPTER 10)');
    CEPslf = CEPrate*hr;    % CEP straight line fit
    plot(hr,CEPslf,'k--','LineWidth',1.5);
else
    disp(['GRAMIAN             = ',num2str(G),' Hr^2.']);
    disp('NO SOLUTION POSSIBLE.');
end;
