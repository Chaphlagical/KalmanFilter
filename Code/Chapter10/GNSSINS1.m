%%  
%%  Grewal and Andrews,
%%  Kalman Filtering: Theory and Practice Using MATLAB,4th Edition,
%%  Wiley, 2014.
%
% COVARIANCE ANALYSIS OF DUAL-FREQUENCY GPS RECEIVER INTEGRATED WITH A
% 1 N.MI./HR INS NAVIGATING ON A 100-KM FIGURE-8 TEST TRACK,
% WITH SIGNAL OUTAGES 1 MINUTE OUT OF EVERY HOUR.
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
EB(1)   = 1e-4;
EB(2)   = 1e-2/3600/180*pi;      % .01 deg/hr/sqrt(s)
EB(3)   = 600;
EB(4)   = 9.8e-5;
EB(5)   = 1e-5;
EB(6)   = 0.03/3600/180*pi;       % 0.03 deg/hr
EB(7)   = 1e-4;
%
%
disp('GPS/INS INTEGRATION');
disp('WITH AN INS HAVING A PARTICULAR ERROR BUDGET');
disp('RESULTING IN 1-NMI/HR CEP RATE');
disp('AND A DUAL-FREQUENCY GPS RECEIVER.');
disp(' ');
disp('IN GENERAL, RESULTS WILL DEPEND ON THE INS ERROR BUDGET');
disp('ALLOCATION DISTRIBUTION, NOT JUST THE RESULTING CEP RATE.');
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
RMSAccNoise     = EB(1);  % [m/s/s/sqrt(s)]
RMSGyroNoise    = EB(2);  % [rad/s/sqrt(s)]
TauSC           = EB(3);  % Correlation time [s] (10 minutes)
RMSAccBias      = EB(4);  % RMS accelerometer bias [m/s/s] (100 micro-G)
RMSAccSF        = EB(5);  % RMS accelerometer scale factor(100 ppm)
RMSGyroBias     = EB(6);  % [rad/s]
RMSGyroSF       = EB(7);  % RMS gyro scale factor(100 ppm)
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
FSC         = -1/TauSC*eye(12);
%
% DYNAMIC DISTURBANCE COVARIANCE MATRIX BLOCKS (ALL TIME-INVARIANT)
%
% Navigation state variables dynamic disturbance noise covariance
%
Qnav        = [zeros(3,9);
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
Q       = [Qnav,zeros(9,14);zeros(12,9),QSC,zeros(12,2);zeros(2,21),Qclock];
%
% INITIAL VALUE OF COVARIANCE MATRIX OF STATE UNCERTAINTY
%
P           = [zeros(9,23);
               zeros(3,9),RMSAccBias^2*eye(3),zeros(3,11);
               zeros(3,12),RMSAccSF^2*eye(3),zeros(3,8);
               zeros(3,15),RMSGyroBias^2*eye(3),zeros(3,5);
               zeros(3,18),RMSGyroSF^2*eye(3),zeros(3,2);
               zeros(2,21),P0clock];
P(1,1)      = 5^2;      % Approx. GNSS-only steady-state value
P(2,2)      = 5^2;      % (ditto)
P(3,3)      = 2^2;      % Altimeter MS error
P(4,4)      = .1^2;     % 0.1 m/s RMS
P(5,5)      = .1^2;     % 0.1 m/s RMS
P(6,6)      = .1^2;     % 0.1 m/s RMS
P(7,7)      = .00001^2; % .01 milliradian RMS
P(8,8)      = .00001^2; % .01 milliradian RMS
P(9,9)      = .00001^2; % .01 milliradian RMS
%
Phi         = [zeros(21,23);zeros(2,21),Phiclock];
%
%
% load 'P0forCEPTest'   % to initialize at low-P0 for CEP tests
%
% Initialize plotting variable arrays
%
CEPINS      = []; % CEP from Riccati equation solution
RMSHorGI    = []; % RMS altitude error with altimeter damping
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
    beforetunnel = (mod(t,3600) < (3600 - 120));
    aftertunnel  = (mod(t,3600) > (3600 - 60));
    if beforetunnel|aftertunnel,
        H           = [HGNSS,zeros(NumSats,18),ones(NumSats,1),zeros(NumSats,1)];
    else
        H           = zeros(NumSats,23);
    end;
    RMSHorGI        = [RMSHorGI,sqrt(P(1,1)+P(2,2))];
    PHT             = P*H';
    K               = PHT/(H*PHT+R);
    P               = P - K*PHT';
    P               = .5*(P + P');
    time            = [time,t];
    RMSHorGI        = [RMSHorGI,sqrt(P(1,1)+P(2,2))];
    [xENU,vENU,aENU,aSensedENU,omegaENU,omegaSensedENU] = Big8TrackSimENU(t,LatRad,LonRad,Alt);
    F11             = Fcore4(LatRad,vENU,aENU);
    F12             = FNSmat12to9(eye(3),aSensedENU,omegaSensedENU);
    PhiINS          = expm(Dt*[F11,F12;zeros(12,9),FSC]);
    Phi(1:21,1:21)  = PhiINS;
    P               = Phi*P*Phi' + Q;
    P               = .5*(P + P');
end;
save 'RMSHorGPSINS' 'RMSHorGI'
save 'TimeSec' 'time'
%
% PLOT BASELINE CEP AND RMS ALTITUDE UNCERTAINTY w/ & \wo ALTIMETER
%
figure;
plot(time/3600,RMSHorGI,'k-','LineWidth',2);
set(gca,'FontSize',14);
%text(1.9,1.9,['CEP RATE = ',num2str(CEPrate),' NMI/HR.'],'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14);
xlabel('TIME [HR]','FontWeight','bold','FontSize',14);
ylabel('HORIZONTAL POSITION UNCERTAINTY [m]','FontWeight','bold','FontSize',14);
title('GPS WITH ~1 NMi/Hr ''MEDIUM ACCURACY'' INS');
%
