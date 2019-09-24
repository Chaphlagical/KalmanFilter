%%  
%%  Grewal and Andrews,
%%  Kalman Filtering: Theory and Practice Using MATLAB,4th Edition,
%%  Wiley, 2014.
%%  
%
% SIMULATION AND COVARIANCE ANALYSIS OF A MEDIUM-ACCURACY INS 
% USING A BAROMETRIC ALTIMETER FOR ALTITUDE STABILIZATION
% WITH SOFTWARE SIMULATED OPERATION ON A 100-KM FIGURE-8 TEST TRACK.
%
close all;
clear all;
%
% Earth model
%
REarth      = 0.6371009e7;  % [m]         mean Earth radius
OmegaEarth  = 0.7292115e-4; % [rad/s]     Earth rotation rate
GM          = 0.3986004e15; % [m^3/s^2]   Earth gravitational constant
%
% Test site location
%
LatDeg      = 39.7931;              % Approximate location of
LonDeg      = -86.2389;             % Indianapolis Motor Speedway,
Alt         = 221;                  % Speedway, IN.
LatRad      = LatDeg*pi/180;        % Latitude in radians
LonRad      = LonDeg*pi/189;        % Longitude in radians
%
% 15-degree horizon mask limit for satellite selection
%
s15 = sin(15*pi/180);               % sine of mask angle
%
% BASELINE INS ERROR BUDGET
%
% INS navigation error budget for 9-state navigation error model with 
% 12-state sensor compensation and barometric altimeter aiding
%
% ERROR BUDGET DATA VECTOR
%
EBname(1,:)   = 'RMS ACCELEROMETER NOISE [m/s/sqrt(s)]           ';
EBname(2,:)   = 'RMS GYRO NOISE [rad/s/sqrt(s)]                  ';
EBname(3,:)   = 'SENSOR COMPENSATION ERROR CORRELATION TIME [s]  ';
EBname(4,:)   = 'RMS ACCELEROMETER BIAS ERROR [m/s/s]            ';
EBname(5,:)   = 'RMS ACCELEROMETER SCALE FACTOR ERROR [part/part]';
EBname(6,:)   = 'RMS GYRO BIAS ERROR [rad/s]                     ';
EBname(7,:)   = 'RMS GYRO SCALE FACTOR ERROR [part/part]         ';
EBname(8,:)   = 'INTERSAMPLE INTERVAL [s]                        ';
EBname(9,:)   = 'RMS ALTIMETER ERROR [m]                         ';
EBname(10,:)  = 'ALTIMETER BIAS CORRELATION TIME [s]             ';
EBname(11,:)  = 'RMS ALTIMETER NOISE [m/sqrt(s)]                 ';
%
EB(1)   = 1e-6;
EB(2)   = 1e-4/3600/180*pi;
EB(3)   = 60;
EB(4)   = 9.8e-5;
EB(5)   = 1e-5;
EB(6)   = 0.01/3600/180*pi;
EB(7)   = 1e-4;
EB(8)   = 1;
EB(9)   = 2;
EB(10)  = 3600;
EB(11)  = 0.1;
%
% DISPLAY BASELINE INS ERROR BUDGET
%
disp('BASELINE INS ERROR BUDGET:');
for item=1:11,
    disp([EBname(item,:),' = ',num2str(EB(item))]);
end;
%
% 9 Navigation Error State Variables + altimeter offset error
%   3 position errors
%   3 veocity errors
%   3 attitude errors
%
%  The state transition matrix depends on dynamic conditions,
%  so it must be built inside the simulation loop.
%
%  The dynamic coupling matrix between sensor compensation errors and
%  INS navigation errors also depends on dynamic conditions.
%
% 12 Sensor Compensation Parameters:
%   3 accelerometer bias errors
%   3 acccelerometer scale factor errors
%   3 gyro bias errors
%   3 gyro scale factor errors
%
% BAROMETRIC ALTIMETER uses local altimeter for ambient error correction
%
TauBaro         = 3600;     % Barometric altimeter error correlation time [s]
RMSBaro         = 2;        % Steady-state RMS barometric altimeter error [m]
RMSAltN         = .1;       % RMS altimeter noise
%
% RMS sensor noise
%
RMSAccNoise     = 1e-6;     % [m/s/s/sqrt(s)]
RMSGyroNoisedph = 1e-4;     % [deg/hr/sqrt(s)]
RMSGyroNoise    = RMSGyroNoisedph/3600/180*pi;  % [rad/s/sqrt(s)]
clear RMSGyroNoisedph;
%
% Sensor compensation errors
%
TauSC           = 600;      % Correlation time [s] (10 minutes)
RMSAccBias      = 9.8e-4;   % RMS accelerometer bias [m/s/s] (100 micro-G)
RMSAccSF        = 1e-4;     % RMS accelerometer scale factor(100 ppm)
RMSGyroBiasdph  = 0.1;      % RMS gyro bias [deg/hr] )
RMSGyroBias     = RMSGyroBiasdph/3600/180*pi;   % [rad/s]
clear RMSGyroBiasdph;
RMSGyroSF       = 1e-4;     % RMS gyro scale factor(100 ppm)
%
RMSAccNoise     = EB(1);  % [m/s/s/sqrt(s)]
RMSGyroNoise    = EB(2);  % [rad/s/sqrt(s)]
TauSC           = EB(3);  % Correlation time [s] (10 minutes)
RMSAccBias      = EB(4);  % RMS accelerometer bias [m/s/s] (100 micro-G)
RMSAccSF        = EB(5);  % RMS accelerometer scale factor(100 ppm)
RMSGyroBias     = EB(6);  % [rad/s]
RMSGyroSF       = EB(7);  % RMS gyro scale factor(100 ppm)
Dt              = EB(8);   % Discrete time-step
RMSBaro         = EB(9);   % Steady-state RMS barometric altimeter error [m]
TauBaro         = EB(10);  % Barometric altimeter error correlation time [s]
RMSAltN         = EB(11);  % RMS altimeter noise
%
% ERROR STATE VARIABLES FOR INS WITH BAROMETRIC ALTIMETER AIDING
%
%   1.  East position error                                 
%   2.  North position error
%   3.  Up position (altitude) error
%   4.  East velocity error
%   5.  North velocity error
%   6.  Up velocity error
%   7.  East-axis attitude rotational (tilt) error
%   8.  North-axis attitude rotational (tilt) error
%   9.  Counterclockwise heading error
%  10.  East accelerometer bias compensation error
%  11.  North accelerometer bias compensation error
%  12.  Vertical accelerometer bias compensation error
%  13.  East accelerometer scale factor compensation error
%  14.  North accelerometer scale factor compensation error
%  15,  Vertical accelerometer scale factor compensation error  
%  16.  East gyro bias compensation error
%  17.  North gyro bias compensation error
%  18.  Vertical gyro bias compensation error
%  19.  East gyro scale factor compensation error
%  20.  North gyro scale factor compensation error
%  21.  Vertical gyro scale factor compensation error
%  22.  Barometer altitude error
%
% DISCRETE TIME STEP (INTERSAMPLE INTERVAL)
%
Dt              = 1;
%  
% TIME-INVARIANT DYNAMIC COEF. MATRIX AND STATE TRANSITION MATRIX BLOCKS
%
FSC         = -1/TauSC*eye(12);
PhiBaro     = exp(-Dt/TauBaro)*eye(1);
Phi         = zeros(22);
Phi(22,22)  = PhiBaro;
%
% DYNAMIC DISTURBANCE COVARIANCE MATRIX BLOCKS (ALL TIME-INVARIANT)
%
% Navigation state variables dynamic disturbance noise covariance
%
Qnav    = [zeros(3,9);
           zeros(3),Dt*RMSAccNoise^2*eye(3),zeros(3);
           zeros(3,6),Dt*RMSGyroNoise^2*eye(3)];
%
% Sensor compensation parameter dynamic disturbance noise covariance
%
QSC     = 2^Dt/TauSC*[RMSAccBias^2*eye(3),zeros(3,9)
                     zeros(3),RMSAccSF^2*eye(3),zeros(3,6)
                     zeros(3,6),RMSGyroBias^2*eye(3),zeros(3);
                     zeros(3,9),RMSGyroSF^2*eye(3)];
%
% Barometric altimeter error dynamic disturbance noise covariance
%
Qbaro   = [2*Dt/TauBaro*RMSBaro^2];
%
% Full 22-state Q-matrix
%
Q       = [Qnav,zeros(9,13);zeros(12,9),QSC,zeros(12,1);zeros(1,21),Qbaro];
%
% ALTIMETER MEASUREMENT SENSITIVITY MATRIX AND MEASUREMENT NOISE COVARIANCE
%
Halt    = [0,0,-1,zeros(1,18),1];  % Altimeter H-matrix
Ralt    = RMSAltN^2;               % Altimeter R-matrix
%
%
% INITIAL VALUE OF COVARIANCE MATRIX OF STATE UNCERTAINTY
%
P           = [zeros(9,22);
               zeros(3,9),RMSAccBias^2*eye(3),zeros(3,10);
               zeros(3,12),RMSAccSF^2*eye(3),zeros(3,7);
               zeros(3,15),RMSGyroBias^2*eye(3),zeros(3,4);
               zeros(3,18),RMSGyroSF^2*eye(3),zeros(3,1);
               zeros(1,21),RMSBaro^2];
P(1,1)      = 5^2;      % Approx. GNSS-only steady-state value
P(2,2)      = 5^2;      % (ditto)
P(3,3)      = 2^2;      % Altimeter MS error
P(4,4)      = .1^2;     % 0.1 m/s RMS
P(5,5)      = .1^2;     % 0.1 m/s RMS
P(6,6)      = .1^2;     % Track RMS value is 1e-5;
P(7,7)      = .0004^2;  % .3 milliradian RMS
P(8,8)      = .0004^2;  % .3 milliradian RMS
P(9,9)      = .0004^2;  % .3 milliradian RMS
%
Pn          = P;        % Covariance without altimeter
P0          = P;        % Save to initialize all cases with the same value
%
% Initialize plotting variable arrays
%
CEPINS      = []; % CEP from Riccati equation solution
RMSalt      = []; % RMS altitude error with altimeter damping
RMSalt0     = []; % RMS altitude error with altimeter damping
time        = []; % time in seconds
%
% Simulated test interval
%
hours   = 8;             % Simulated test interval
minutes = 60*hours;
seconds = 60*minutes;
%
for t=0:Dt:seconds,
    time            = [time,t];
    CEP             = 1.2*sqrt(P(1,1)+P(2,2))/1852; % CEP in nautical miles
    CEPINS          = [CEPINS,CEP];
    RMSalt          = [RMSalt,sqrt(P(3,3))];
    RMSalt0         = [RMSalt0,sqrt(Pn(3,3))];
    PHT             = P*Halt';
    K               = PHT/(Halt*PHT+Ralt);
    P               = P - K*PHT';
    P               = .5*(P + P');
    time            = [time,t];
    CEP             = 1.2*sqrt(P(1,1)+P(2,2))/1852; % CEP in nautical miles
    CEPINS          = [CEPINS,CEP];
    RMSalt          = [RMSalt,sqrt(P(3,3))];
    RMSalt0         = [RMSalt0,sqrt(Pn(3,3))];
    [xENU,vENU,aENU,aSensedENU,omegaENU,omegaSensedENU] = Big8TrackSimENU(t,LatRad,LonRad,Alt);
    F11             = Fcore4(LatRad,vENU,aENU);
    F12             = FNSmat12to9(eye(3),aSensedENU,omegaSensedENU);
    PhiINS          = expm(Dt*[F11,F12;zeros(12,9),FSC]);
    Phi(1:21,1:21)  = PhiINS;
    P               = Phi*P*Phi' + Q;
    P               = .5*(P + P');
    Pn              = Phi*Pn*Phi' + Q;
    Pn              = .5*(Pn + Pn');
end;
%
% LEAST-SQUARES FIT CEP RATE
%
hr          = time/3500;
CEPrate     =  hr*CEPINS'/(hr*hr');
disp(['BASELINE CEP RATE = ',num2str(CEPrate),' NMI/HR.']);
%
% PLOT BASELINE CEP AND RMS ALTITUDE UNCERTAINTY w/ & \wo ALTIMETER
%
figure;
plot(hr,CEPINS,'k-','LineWidth',2);
set(gca,'FontSize',14);
%text(1.9,1.9,['CEP RATE = ',num2str(CEPrate),' NMI/HR.'],'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
xlabel('TIME [HR]','FontWeight','bold','FontSize',14);
ylabel('CEP [NMi]','FontWeight','bold','FontSize',14);
title('100 KM FIGURE-8 TRACK, INS PLUS ALTIMETER');
%
figure;
semilogy(hr,RMSalt0,'k--',hr,RMSalt,'k-','LineWidth',2);
set(gca,'FontSize',14);
legend('INERTIAL ONLY','ALTIMETER AIDING');
xlabel('TIME [HR]','FontWeight','bold','FontSize',14);
ylabel('RMS ALTITUDE ERROR [m]','FontWeight','bold','FontSize',14);
title('100 KM FIGURE-8 TRACK, INS \pm ALTIMETER');

