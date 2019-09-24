%%  
%%  Grewal and Andrews,
%%  Kalman Filtering: Theory and Practice Using MATLAB,4th Edition,
%%  Wiley, 2014.
%
% COVARIANCE ANALYSIS OF DUAL-FREQUENCY GPS RECEIVER INTEGRATED WITH A
% 0.1, 10, & 10 N.MI./HR INS NAVIGATING ON A 100-KM FIGURE-8 TEST TRACK,
% WITH SIGNAL OUTAGE FOR LAST HOUR.
%
% STATE VECTOR:
%   9 INS ERROR MODEL DYNAMIC VARIABLES
%  12 INERTIAL SENSOR COMPENSATION PARAMETERS
%   2 RECEIVER CLOCK ERROR VERIABLES
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
%%%%    
%%%%    %%%%%%%%    %%%%    %%%%%        %%%%%%                %%%    
%%%%    %%%%%%%%    %%%%%   %%%%%      %%%%%%%%%%             %%%%    
%%%%      %%%%       %%%%%  %%%%     %%%%%%   %%%%%          %%%%%         
%%%%      %%%%       %%%%%% %%%%     %%%%%                     %%%
%%%%      %%%%       %%%% %%%%%%      %%%%%%%%%%%%             %%%
%%%%      %%%%       %%%%  %%%%%        %%%%%%%%%%%            %%%
%%%%      %%%%       %%%%   %%%%               %%%%            %%%
%%%%    %%%%%%%%    %%%%%    %%%%%    %%%%%   %%%%            %%%%%
%%%%    %%%%%%%%    %%%%%     %%%%      %%%%%%%%%             %%%%%   
%%%%    
%
% INS navigation error budget for 9-state navigation error model with 
% 12-state sensor compensation and barometric altimeter aiding
%
% INS ERROR BUDGET DESCRIPTIONS
%
EBname(1,:)   = 'RMS ACCELEROMETER NOISE [m/s/sqrt(s)]           ';
EBname(2,:)   = 'RMS GYRO NOISE [rad/s/sqrt(s)]                  ';
EBname(3,:)   = 'SENSOR COMPENSATION ERROR CORRELATION TIME [s]  ';
EBname(4,:)   = 'RMS ACCELEROMETER BIAS ERROR [m/s/s]            ';
EBname(5,:)   = 'RMS ACCELEROMETER SCALE FACTOR ERROR [part/part]';
EBname(6,:)   = 'RMS GYRO BIAS ERROR [rad/s]                     ';
EBname(7,:)   = 'RMS GYRO SCALE FACTOR ERROR [part/part]         ';
%
% 10 NMi/Hr Error Budget
%
EB1(1)   = 1e-4;
EB1(2)   = 1e-2/3600/180*pi;      % .01 deg/hr/sqrt(s)
EB1(3)   = 1800;
EB1(4)   = 9.8e-4;
EB1(5)   = 1e-4;
EB1(6)   = 0.19/3600/180*pi;       % 0.19 deg/hr
EB1(7)   = 1e-3;                   % 1e-3 for 10 NMi/Hr CEP rate
%
% 1 NMi/Hr Error Budget
%
EB2(1)   = 1e-4;
EB2(2)   = 1e-2/3600/180*pi;      % .01 deg/hr/sqrt(s)
EB2(3)   = 600;
EB2(4)   = 9.8e-5;
EB2(5)   = 1e-5;
EB2(6)   = 0.03/3600/180*pi;       % 0.03 deg/hr
EB2(7)   = 1e-4;
%
% 0.1 NMi/Hr Error Budget
%
EB3(1)   = 1e-5;
EB3(2)   = 1e-4/3600/180*pi;      % .0001 deg/hr/sqrt(s)
EB3(3)   = 15;                    % 15 for 0.1 NMi/Hr CEP Rate 1800
EB3(4)   = 9.8e-5;
EB3(5)   = 1e-4;
EB3(6)   = 0.01/3600/180*pi;      % 0.01/3600/180*pi for 0.1 NMi/Hr CEP Rate
EB3(7)   = 1e-4;                  % 1e-4 for 1 NMi/Hr CEP rate
%
disp('GPS/INS INTEGRATION');
disp('WITH LOW-, MEDIUM-, AND HIGH-ACCURACY INSs:');
disp('10 N.Mi/Hr, 1 N.Mi/Hr, 0.1 N.Mi/Hr CEP RATES');
disp('INTEGRATED WITH DUAL-FREQUENCY GPS RECEIVER.');
disp(' ');
disp('(IN GENERAL, RESULTS WILL DEPEND ON THE INS ERROR BUDGET');
disp('ALLOCATION DISTRIBUTION, NOT JUST THE RESULTING CEP RATE.)');
%
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
% INITIALIZE ERROR BUDGET VARIABLES
%
RMSAccNoise1     = EB1(1);  % [m/s/s/sqrt(s)]
RMSGyroNoise1    = EB1(2);  % [rad/s/sqrt(s)]
TauSC1           = EB1(3);  % Correlation time [s] (10 minutes)
RMSAccBias1      = EB1(4);  % RMS accelerometer bias [m/s/s] (100 micro-G)
RMSAccSF1        = EB1(5);  % RMS accelerometer scale factor(100 ppm)
RMSGyroBias1     = EB1(6);  % [rad/s]
RMSGyroSF1       = EB1(7);  % RMS gyro scale factor(100 ppm)
%
RMSAccNoise2     = EB2(1);  % [m/s/s/sqrt(s)]
RMSGyroNoise2    = EB2(2);  % [rad/s/sqrt(s)]
TauSC2           = EB2(3);  % Correlation time [s] (10 minutes)
RMSAccBias2      = EB2(4);  % RMS accelerometer bias [m/s/s] (100 micro-G)
RMSAccSF2        = EB2(5);  % RMS accelerometer scale factor(100 ppm)
RMSGyroBias2     = EB2(6);  % [rad/s]
RMSGyroSF2       = EB2(7);  % RMS gyro scale factor(100 ppm)
%
RMSAccNoise3     = EB3(1);  % [m/s/s/sqrt(s)]
RMSGyroNoise3    = EB3(2);  % [rad/s/sqrt(s)]
TauSC3           = EB3(3);  % Correlation time [s] (10 minutes)
RMSAccBias3      = EB3(4);  % RMS accelerometer bias [m/s/s] (100 micro-G)
RMSAccSF3        = EB3(5);  % RMS accelerometer scale factor(100 ppm)
RMSGyroBias3     = EB3(6);  % [rad/s]
RMSGyroSF3       = EB3(7);  % RMS gyro scale factor(100 ppm)
Dt              =1;       % Discrete time-step
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
%  22.  Receiver clock bias [m]
%  23.  Receiver clock drift [m/s]
%   
% TIME-INVARIANT DYNAMIC COEF. SUB-MATRIX
%
FSC1        = -1/TauSC1*eye(12);
FSC2        = -1/TauSC2*eye(12);
FSC3        = -1/TauSC3*eye(12);
%
% DYNAMIC DISTURBANCE COVARIANCE MATRIX BLOCKS (ALL TIME-INVARIANT)
%
% Navigation state variables dynamic disturbance noise covariance
%
Qnav1        = [zeros(3,9);
           zeros(3),Dt*RMSAccNoise1^2*eye(3),zeros(3);
           zeros(3,6),Dt*RMSGyroNoise1^2*eye(3)];
%
Qnav2        = [zeros(3,9);
           zeros(3),Dt*RMSAccNoise2^2*eye(3),zeros(3);
           zeros(3,6),Dt*RMSGyroNoise2^2*eye(3)];
%
Qnav3        = [zeros(3,9);
           zeros(3),Dt*RMSAccNoise3^2*eye(3),zeros(3);
           zeros(3,6),Dt*RMSGyroNoise3^2*eye(3)];
%
% Sensor compensation parameter dynamic disturbance noise covariance
%
QSC1     = 2^Dt/TauSC1*[RMSAccBias1^2*eye(3),zeros(3,9)
                     zeros(3),RMSAccSF1^2*eye(3),zeros(3,6)
                     zeros(3,6),RMSGyroBias1^2*eye(3),zeros(3);
                     zeros(3,9),RMSGyroSF1^2*eye(3)];
%
QSC2     = 2^Dt/TauSC2*[RMSAccBias2^2*eye(3),zeros(3,9)
                     zeros(3),RMSAccSF2^2*eye(3),zeros(3,6)
                     zeros(3,6),RMSGyroBias2^2*eye(3),zeros(3);
                     zeros(3,9),RMSGyroSF2^2*eye(3)];
%
QSC3     = 2^Dt/TauSC3*[RMSAccBias3^2*eye(3),zeros(3,9)
                     zeros(3),RMSAccSF3^2*eye(3),zeros(3,6)
                     zeros(3,6),RMSGyroBias3^2*eye(3),zeros(3);
                     zeros(3,9),RMSGyroSF3^2*eye(3)];
%
%
% Set simulation time-step
%
Dt          = 1;  % 1-second discrete time-step
%
%%%%   
%%%%         %%%%%%%        %%%%%%%%%%           %%%%%%
%%%%      %%%%%%%%%%%%      %%%%%%%%%%%%%      %%%%%%%%%%
%%%%    %%%%%     %%%%%%     %%%%   %%%%%     %%%%%   %%%%%
%%%%   %%%%                  %%%%   %%%%%     %%%%%
%%%%   %%%%     %%%%%%%%     %%%%%%%%%%%      %%%%%%%%%%%%
%%%%   %%%%     %%%%%%%%     %%%%%%%%%          %%%%%%%%%%%
%%%%    %%%%%      %%%       %%%%                      %%%%
%%%%      %%%%%%%%%%%%%%    %%%%%%            %%%%%   %%%%
%%%%         %%%%%%    %%   %%%%%%              %%%%%%%%%
%%%%
%%%%  GPS SINGLE-FREQUENCY RECEIVER MODEL
%%%%
%%%%
%  
% Set up the GNSS (GPS) simulation database
%
% SINGLE-FREQUENCY RECEIVER REQUIRES 2-STATE CLOCK MODEL
% PLUS ESPONENTIALLY CORRELATED IONO DELAY ERROR STATE FOR EACH SATELLITE
%
% RECEIVER CLOCK MODEL
%
% Receiver clock model parameters
%
taudrift    = 3600;     % Clock drift rate correlation time [1 hr]
sigmadrift  = 10;       % RMS clock drift rate [m/s]
[Fclock,Qclock,Hclock,P0clock] = ClockModel2(sigmadrift,taudrift);
Phiclock    = expm(Fclock*Dt); % state transition matrix
Qclock      = Dt*Qclock;
sigmaclock  = sqrt(Qclock(2,2));
%
% GPS SATELLITE EPHEMERIDES (USED FOR SIMULATION)
%
% FOR SIMPLICITY IN COVARIANCE ANALYSIS, GPS SATELLITE SIMULATION
% APPROXIMATES ALL SATELLITE ORBITS AS CIRCULAR WITH 55 DEG INCLINATION.
% THE POSITION IN ORBIT IS THEN CHARACTERIZED BY RIGHT ASCENSION OF THE
% ASCENDING NODE (RA) AND THE IN-PLANE ARGUMENT (PA) OF THE SATELLITE
% MEASURED FROM THE ASCENDING NODE AT TIME t=0.
% 
% YUMAdata loads right ascension (RA) and perigee argument (PA) only.
% The satellite simulator HSatENU.m assumes all orbits are circular.
%
global RA PA ;% globalize satellite ephemerides for use by functions
YUMAdata;
%
% Check for data mismatch
%
if length(RA) ~= length(PA), error('Ephemeris data corrupted!'); return;end;
% 
%
% GPS SATELLITE SIGNAL ERROR MODEL PARAMETERS FOR SINGLE-FREQUENCY RECEIVER
%
NumSats         = length(RA);  % number of satellites listed in ephemerides
%
% RECEIVER PSEUDORANGE MEASUREMENT NOISE
%
sigmaPRMN       = 10;  % RMS pseudorange measurement noise [m]
%
% SATELLITE SIGNAL PSEUDORANGE NOISE COVARIANCE MATRIX
%
R               = sigmaPRMN^2*eye(NumSats);
%
%
% Full 22-state Q-matrix
%
Q1       = [Qnav1,zeros(9,14);zeros(12,9),QSC1,zeros(12,2);zeros(2,21),Qclock];
%
Q2       = [Qnav2,zeros(9,14);zeros(12,9),QSC2,zeros(12,2);zeros(2,21),Qclock];
%
Q3       = [Qnav3,zeros(9,14);zeros(12,9),QSC3,zeros(12,2);zeros(2,21),Qclock];
%
% INITIAL VALUE OF COVARIANCE MATRIX OF STATE UNCERTAINTY
%
P1           = [zeros(9,23);
               zeros(3,9),RMSAccBias1^2*eye(3),zeros(3,11);
               zeros(3,12),RMSAccSF1^2*eye(3),zeros(3,8);
               zeros(3,15),RMSGyroBias1^2*eye(3),zeros(3,5);
               zeros(3,18),RMSGyroSF1^2*eye(3),zeros(3,2);
               zeros(2,21),P0clock];
P1(1,1)      = 2^2;      % Approx. GNSS-only steady-state value
P1(2,2)      = 2^2;      % (ditto)
P1(3,3)      = 2^2;      % Altimeter MS error
P1(4,4)      = .1^2;     % 0.1 m/s RMS
P1(5,5)      = .1^2;     % 0.1 m/s RMS
P1(6,6)      = .1^2;     % 0.1 m/s RMS
P1(7,7)      = .00001^2; % .01 milliradian RMS
P1(8,8)      = .00001^2; % .01 milliradian RMS
P1(9,9)      = .00001^2; % .01 milliradian RMS
%
P2           = [zeros(9,23);
               zeros(3,9),RMSAccBias2^2*eye(3),zeros(3,11);
               zeros(3,12),RMSAccSF2^2*eye(3),zeros(3,8);
               zeros(3,15),RMSGyroBias2^2*eye(3),zeros(3,5);
               zeros(3,18),RMSGyroSF2^2*eye(3),zeros(3,2);
               zeros(2,21),P0clock];
P2(1,1)      = 2^2;      % Approx. GNSS-only steady-state value
P2(2,2)      = 2^2;      % (ditto)
P2(3,3)      = 2^2;      % Altimeter MS error
P2(4,4)      = .1^2;     % 0.1 m/s RMS
P2(5,5)      = .1^2;     % 0.1 m/s RMS
P2(6,6)      = .1^2;     % 0.1 m/s RMS
P2(7,7)      = .00001^2; % .01 milliradian RMS
P2(8,8)      = .00001^2; % .01 milliradian RMS
P2(9,9)      = .00001^2; % .01 milliradian RMS
%
P3           = [zeros(9,23);
               zeros(3,9),RMSAccBias3^2*eye(3),zeros(3,11);
               zeros(3,12),RMSAccSF3^2*eye(3),zeros(3,8);
               zeros(3,15),RMSGyroBias3^2*eye(3),zeros(3,5);
               zeros(3,18),RMSGyroSF3^2*eye(3),zeros(3,2);
               zeros(2,21),P0clock];
P3(1,1)      = 2^2;      % Approx. GNSS-only steady-state value
P3(2,2)      = 2^2;      % (ditto)
P3(3,3)      = 2^2;      % Altimeter MS error
P3(4,4)      = .1^2;     % 0.1 m/s RMS
P3(5,5)      = .1^2;     % 0.1 m/s RMS
P3(6,6)      = .1^2;     % 0.1 m/s RMS
P3(7,7)      = .00001^2; % .01 milliradian RMS
P3(8,8)      = .00001^2; % .01 milliradian RMS
P3(9,9)      = .00001^2; % .01 milliradian RMS
%
Phi1        = [zeros(21,23);zeros(2,21),Phiclock];
Phi2        = [zeros(21,23);zeros(2,21),Phiclock];
Phi3        = [zeros(21,23);zeros(2,21),Phiclock];
%
%
% load 'P0forCEPTest'   % to initialize at low-P0 for CEP tests
%
% Initialize plotting variable arrays
%
RMSHorGI1   = []; % RMS horizontal error, 10 NMI/Hr
RMSHorGI2   = []; % RMS horizontal error, 1 NMi/Hr
RMSHorGI3   = []; % RMS horizontal error, 0.1 NMi/Hr
time        = []; % time in seconds
%
% Simulated test interval
%
hours   = 4;             % Simulated test interval
minutes = 60*hours;
seconds = 60*minutes;
%
for t=0:Dt:seconds,
    time            = [time,t];
    [xENU,vENU,aENU,aSensedENU,omegaENU,omegaSensedENU] = Big8TrackSimENU(t,LatRad,LonRad,Alt); 
    HGNSS   = HSatENU(t,LatRad+xENU(2)/REarth,LonRad,Alt);
    NoSatsAvail = NumSats;
    for j=1:NumSats,
        if HGNSS(j,3) > -0.25881904510252 
            NoSatsAvail = NoSatsAvail - 1;
            HGNSS(j,1) = 0;
            HGNSS(j,2) = 0;
            HGNSS(j,3) = 0;
        end;
    end;
    RMSHorGI1       = [RMSHorGI1,sqrt(P1(1,1)+P1(2,2))];
    RMSHorGI2       = [RMSHorGI2,sqrt(P2(1,1)+P2(2,2))];
    RMSHorGI3       = [RMSHorGI3,sqrt(P3(1,1)+P3(2,2))];
    if t<= seconds - 3600,
        H               = [HGNSS,zeros(NumSats,18),ones(NumSats,1),zeros(NumSats,1)];
        PHT1            = P1*H';
        PHT2            = P2*H';
        PHT3            = P3*H';
        K1              = PHT1/(H*PHT1+R);
        K2              = PHT2/(H*PHT2+R);
        K3              = PHT3/(H*PHT3+R);
        P1              = P1 - K1*PHT1';
        P1              = .5*(P1 + P1');
        P2              = P2 - K2*PHT2';
        P2              = .5*(P2 + P2');
        P3              = P3 - K3*PHT3';
        P3              = .5*(P3 + P3');
    end;
    time            = [time,t];
    RMSHorGI1       = [RMSHorGI1,sqrt(P1(1,1)+P1(2,2))];
    RMSHorGI2       = [RMSHorGI2,sqrt(P2(1,1)+P2(2,2))];
    RMSHorGI3       = [RMSHorGI3,sqrt(P3(1,1)+P3(2,2))];
    [xENU,vENU,aENU,aSensedENU,omegaENU,omegaSensedENU] = Big8TrackSimENU(t,LatRad,LonRad,Alt);
    F11             = Fcore4(LatRad,vENU,aENU);
    F12             = FNSmat12to9(eye(3),aSensedENU,omegaSensedENU);
    PhiINS1         = expm(Dt*[F11,F12;zeros(12,9),FSC1]);
    Phi1(1:21,1:21) = PhiINS1;
    PhiINS2         = expm(Dt*[F11,F12;zeros(12,9),FSC2]);
    Phi2(1:21,1:21) = PhiINS2;
    PhiINS3         = expm(Dt*[F11,F12;zeros(12,9),FSC3]);
    Phi3(1:21,1:21) = PhiINS3;
    P1              = Phi1*P1*Phi1' + Q1;
    P1              = .5*(P1 + P1');
    P2              = Phi2*P2*Phi2' + Q2;
    P2              = .5*(P2 + P2');
    P3              = Phi3*P3*Phi3' + Q3;
    P3              = .5*(P3 + P3');
end;
%save 'RMSHorGPSINS' 'RMSHorGI'
%save 'TimeSec' 'time'
%
% PLOT BASELINE RMS HORIZONTAL UNCERTAINTIES FOR ALL THREE MODELS
%
figure;
semilogy(time/3600,RMSHorGI1,'k-',time/3600,RMSHorGI2,'k--',time/3600,RMSHorGI3,'k:','LineWidth',2);
set(gca,'FontSize',14);
%text(1.9,1.9,['CEP RATE = ',num2str(CEPrate),' NMI/HR.'],'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
xlabel('TIME [HR]','FontWeight','bold','FontSize',14);
ylabel('RMS HORIZONTAL POS. UNCERTAINTY [m]','FontWeight','bold','FontSize',14);
title('DUAL-FREQ. GPS WITH 10, 1, AND 0.1 NMi/Hr INS');
%
