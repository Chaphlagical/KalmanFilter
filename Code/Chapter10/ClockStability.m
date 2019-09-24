%
% Grewal & Andrews, Kalman Filtering: Theory and Practice Using MATLAB,
% 4th Edition, Wiley, 2014
%
%
% Influence of receiver clock stability on timing accuracy
%
% Varies RMS clock stability over 1-second interval from 10^(-9) to 10^(-8)
% Simulates 100 minutes of filtering with that geometry.
% Plots RMS position and clock uncertainties versus clock stability
%
% Assumes 15 meter RMS pseudorange error, "good" satellite geometry
%
close all,
clear all,
%
Dt          = 1;  % Discrete time-step
%
P0          = zeros(5);
P0(1,1)     = 1e6;
P0(2,2)     = 1e6;
P0(3,3)     = 1e6;
P0(4,4)     = 9e6;
P0(5,5)     = 9e2; % Initial covariance of state estimation uncertainty
%
Phi = zeros(5);
Phi(1,1)    = 1;
Phi(2,2)    = 1;
Phi(3,3)    = 1;
Phi(4,4)    = 1;
Phi(4,5)    = Dt;
Phi(5,5)    = 1; % State transition matrix
%
% Dynamic disturbance noise covariance for sigma_f/f = 1
%
Q1          = zeros(5);
Q1(4,4)     = 6e16;
Q1(4,5)     = 5e16;
Q1(5,4)     = Q1(4,5);
Q1(5,5)     = 10e16; 
%
R           = 225*eye(4); % Covariance of measurement noise
r           = 225;        % Diagonal term for serial processing
%
% Four satellites with good geometry
%
Elev(1)     = 0;
Azim(1)     = 0;
Elev(2)     = 0;
Azim(2)     = 120;
Elev(3)     = 0;
Azim(3)     = 240;
Elev(4)     = 90;
Azim(4)     = 0;
H           = zeros(4,5);
H(1,4)      = 1;
H(2,4)      = 1;
H(3,4)      = 1;
H(4,4)      = 1;
phi1        = Elev(1)*pi/180;
sphi1       = sin(phi1);
cphi1       = cos(phi1);
H(1,3)      = sphi1;
theta1      = Azim(1)*pi/180;
stheta1     = sin(theta1);
ctheta1     = cos(theta1);
H(1,1)      = -ctheta1*cphi1;
H(1,2)      = -stheta1*cphi1;
phi2        = Elev(2)*pi/180;
sphi2       = sin(phi2);
cphi2       = cos(phi2);
H(2,3)      = sphi2;
theta2      = Azim(2)*pi/180;
stheta2     = sin(theta2);
ctheta2     = cos(theta2);
H(2,1)      = -ctheta2*cphi2;
H(2,2)      = -stheta2*cphi2;
phi3        = Elev(3)*pi/180;
sphi3       = sin(phi3);
cphi3       = cos(phi3);
H(3,3)      = sphi3;
theta3      = Azim(3)*pi/180;
stheta3     = sin(theta3);
ctheta3     = cos(theta3);
H(3,1)      = -ctheta3*cphi3;
H(3,2)      = -stheta3*cphi3;
phi4        = Elev(4)*pi/180;
sphi4       = sin(phi4);
cphi4       = cos(phi4);
H(4,3)      = sphi4;
theta4      = Azim(4)*pi/180;
stheta4     = sin(theta4);
ctheta4     = cos(theta4);
H(4,1)      = -ctheta4*cphi4;
H(4,2)      = -stheta4*cphi4;
for k=1:4,
    RMSclockSt(k)  = 10^(k-10);
    Q           = RMSclockSt(k)^2*Q1;    
    P        = P0;
    for j=1:4,
        h   = H(j,:);
        PHT = P*h';
        K   = PHT/(h*PHT+r);
        P   = P - K*PHT';
    end;
    P   = (P+P')/2;
    %
    for m=1:3600, % ten minutes of filtering
        P   = Phi*P*Phi'+Q;
        for j=1:4,
            h   = H(j,:);
            PHT = P*h';
            K   = PHT/(h*PHT+r);
            P   = P - K*PHT';
        end;
        P   = (P+P')/2;
    end;
    NoErr(k) = sqrt(P(1,1));
    EaErr(k) = sqrt(P(2,2));
    DoErr(k) = sqrt(P(3,3));
    CbErr(k) = sqrt(P(4,4));
    CdErr(k) = sqrt(P(5,5));
end;
figure;
semilogx(RMSclockSt,NoErr,'k-',RMSclockSt,EaErr,'k*',RMSclockSt,DoErr,'k:');
title('RMS Navigation Uncertainty vs Clock Stability');
xlabel('Clock Stability [1-sec, part/part]')
ylabel('RMS after 1 Hour [m]');
legend('North','East','Altitude');
%
figure;
loglog(RMSclockSt,CbErr,'k-',RMSclockSt,CdErr,'k--','LineWidth',2);
title('Clock Uncertainty vs Clock Stability','FontSize',14);
xlabel('RMS Clock Stability [part/part @ 1 s]','FontSize',14)
ylabel('RMS Uncertainty after 1 Hour of Filtering','FontSize',14);
legend('C_b [m]','C_d [m/s]');
