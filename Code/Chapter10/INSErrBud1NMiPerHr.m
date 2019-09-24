%%  
%%  Grewal and Andrews,
%%  Kalman Filtering: Theory and Practice Using MATLAB,4th Edition,
%%  Wiley, 2014.
%
% PLEASE READ CHAPTER 10 FOR CAVEATS RELATED TO ERROR BUDGET M-FILES
%  
% GENERATES AN 11-PARAMETER ERROR BUDGET FOR A MEDUIM-ACCURACY INERTIAL
% NAVIGATION SYSTEM WITH ABOUT 1 NAUTICAL MILE PER HOUR CEP RATE, USING 
% A BAROMETRIC ALTIMETER FOR VERTICAL CHANNEL (ALTITUDE) STABILIZATION.
%
% PERFORMANCE IS FIRST VERIFIED BY SIMULATED OPERATION ON A 100-KM
% FIGURE-8 TEST TRACK FOR 4 HOURS.  CEP IS PLOTTED IN NAUTICAL MILES
% VERSUS TIME IN HOURS, LEAST-SQUARES-FIT CEP IS DISPLAYED, AND RMS 
% ALTITUDE ERRROR IS PLOTTED IN METERS---WITH AND WITHOUT ALTIMETER AIDING 
%
% ERROR BUDGET CAN BE saveD TO FILE 'INSEBMedAcc.mat'
%
% NEXT PERFORMS DIFFERENTIAL  ANALYSIS OF THE DEPENDENCE OF CEP RATE ON
% EACH EACH OF THE 9 ERROR BUDGET PARAMETERS SPECIFIC TO THE INS.
%
% CEP RATE IS ESTIMATED FOR EACH OF THE 9 PARAMETERS SCALED UP AND DOWN BY
% AN ORDER OF MAGNITUDE, ESTIMATING CEP RATE BY SIMULATED OPERATION ON A 
% 100-KM FIGURE-8 TEST TRACK.
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
EB(1)   = 1e-4;
EB(2)   = 1e-2/3600/180*pi;      % .01 deg/hr/sqrt(s)
EB(3)   = 600;
EB(4)   = 9.8e-5;
EB(5)   = 1e-5;
EB(6)   = 0.03/3600/180*pi;       % 0.03/3600/180*pi for 1 NMi/Hr CEP Rate
EB(7)   = 1e-4;                   % 1e-4 for 1 NMi/Hr CEP rate
EB(8)   = 1;
EB(9)   = 2;
EB(10)  = 3600;
EB(11)  = 0.1;
%
%
% DISPLAY BASELINE INS ERROR BUDGET
%
disp('BASELINE INS ERROR BUDGET:');
for item=1:11,
    disp([EBname(item,:),' = ',num2str(EB(item))]);
end;
disp(['RMS GYRO NOISE [deg/hr/sqrt(s)]                  = ',num2str(EB(2)*3600*180/pi)]);
disp(['RMS GYRO BIAS ERROR [deg/hr]                     = ',num2str(EB(6)*3600*180/pi)]);
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
P(7,7)      = .00001^2; % .4 milliradian RMS
P(8,8)      = .00001^2; % .4 milliradian RMS
P(9,9)      = .00001^2; % .4 milliradian RMS
%
Pn          = P;        % Covariance without altimeter
P0          = P;        % Save to initialize all cases with the same value
%
%
% load 'P0forCEPTest'   % to initialize at low-P0 for CEP tests
%
% Initialize plotting variable arrays
%
CEPINS      = []; % CEP from Riccati equation solution
RMSHorINS   = [];
RMSalt      = []; % RMS altitude error with altimeter damping
RMSalt0     = []; % RMS altitude error with altimeter damping
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
    CEP             = 1.2*sqrt(P(1,1)+P(2,2))/1852; % CEP in nautical miles
    CEPINS          = [CEPINS,CEP];
    RMSHorINS       = [RMSHorINS,sqrt(P(1,1)+P(2,2))];
    RMSalt          = [RMSalt,sqrt(P(3,3))];
    RMSalt0         = [RMSalt0,sqrt(Pn(3,3))];
    PHT             = P*Halt';
    K               = PHT/(Halt*PHT+Ralt);
    P               = P - K*PHT';
    P               = .5*(P + P');
    time            = [time,t];
    CEP             = 1.2*sqrt(P(1,1)+P(2,2))/1852; % CEP in nautical miles
    CEPINS          = [CEPINS,CEP];
    RMSHorINS       = [RMSHorINS,sqrt(P(1,1)+P(2,2))];
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
save 'RMSHorINS1m' 'RMSHorINS'
hr          = time/3600;
figure;
plot(hr,RMSHorINS,'k-','LineWidth',2);
set(gca,'FontSize',14);
xlabel('TIME [HR]','FontWeight','bold','FontSize',14);
ylabel('HORIZONTAL POSITION UNCERTAINTY [m]','FontWeight','bold','FontSize',14);
title('~1 NMi/Hr ''MEDIUM ACCURACY'' INS');
%
% LEAST-SQUARES FIT CEP RATE
%
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
title('~1 NMi/Hr ''MEDIUM ACCURACY'' INS');
%
figure;
semilogy(hr,RMSalt0,'k--',hr,RMSalt,'k-','LineWidth',2);
set(gca,'FontSize',14);
legend('INERTIAL ONLY','ALTIMETER AIDING');
xlabel('TIME [HR]','FontWeight','bold','FontSize',14);
ylabel('RMS ALTITUDE ERROR [m]','FontWeight','bold','FontSize',14);
title('~1 NMi/Hr ''MEDIUM ACCURACY'' INS');
%
%
%
response = input('Save error budget to file ''INSEBMedAcc.mat'' ?  (Yes or No): ','s');
if (response(1)=='Y')|(response(1)=='y')
    save 'INSEBMedAcc' 'EB';
    disp('Error budget saved to file INSEBMedAcc.mat');
end;
%
% COMPUTE THE CEP FOR 1st 7 INS BUDGET PARAMETERS SCALED UP AND DOWN BY 10
%
CEPrate0    = CEPrate;  % BASELINE CEEP
CEPrate1    = [];       % BASELINE w PARAMETER SCALED DOWN BY 10
CEPrate2    = [];       % BASELINE w PARAMETER SCALED UP BY 10
EBindex     = (1:7);    % INDICES OF FIRST 7 PARAMETERS
row = 0;
for LogScale=-1:2:1;
    Scale   = 10^LogScale;
    row     = row + 1;
    for BudgetItem=1:7,
        SEB             = EB;
        SEB(BudgetItem) = Scale*EB(BudgetItem);
        RMSAccNoise     = SEB(1);  % [m/s/s/sqrt(s)]
        RMSGyroNoise    = SEB(2);  % [rad/s/sqrt(s)]
        TauSC           = SEB(3);  % Correlation time [s] (10 minutes)
        RMSAccBias      = SEB(4);  % RMS accelerometer bias [m/s/s] (100 micro-G)
        RMSAccSF        = SEB(5);  % RMS accelerometer scale factor(100 ppm)
        RMSGyroBias     = SEB(6);  % [rad/s]
        RMSGyroSF       = SEB(7);  % RMS gyro scale factor(100 ppm)
        Dt              = EB(8);   % Discrete time-step
        RMSBaro         = EB(9);   % Steady-state RMS barometric altimeter error [m]
        TauBaro         = EB(10);  % Barometric altimeter error correlation time [s]
        RMSAltN         = EB(11);  % RMS altimeter noise
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
        P           = P0;
        %
        %
        % Initialize plotting variable arrays
        %
        CEPINS      = []; % CEP from Riccati equation solution
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
            CEP             = 1.2*sqrt(P(1,1)+P(2,2))/1852; % CEP in nautical miles
            CEPINS          = [CEPINS,CEP];
            PHT             = P*Halt';
            K               = PHT/(Halt*PHT+Ralt);
            P               = P - K*PHT';
            P               = .5*(P + P');
            time            = [time,t];
            CEP             = 1.2*sqrt(P(1,1)+P(2,2))/1852; % CEP in nautical miles
            CEPINS          = [CEPINS,CEP];
            [xENU,vENU,aENU,aSensedENU,omegaENU,omegaSensedENU] = Big8TrackSimENU(t,LatRad,LonRad,Alt);
            F11             = Fcore4(LatRad,vENU,aENU);
            F12             = FNSmat12to9(eye(3),aSensedENU,omegaSensedENU);
            PhiINS          = expm(Dt*[F11,F12;zeros(12,9),FSC]);
            Phi(1:21,1:21)  = PhiINS;
            P               = Phi*P*Phi' + Q;
            P               = .5*(P + P');
        end;
        hr          = time/3600;
        CEPrate     =  hr*CEPINS'/(hr*hr');
        if Scale==0.1
            CEPrate1    = [CEPrate1,CEPrate];
        else
            CEPrate2    = [CEPrate2,CEPrate];
        end;
    end;
end;
%
% PLOT RESULTS
%
figure;
semilogy(EBindex,CEPrate2/CEPrate0,'k^',EBindex,CEPrate1/CEPrate0,'kv',.5,.1,'w.','LineWidth',2);
set(gca,'FontSize',14);
legend('SCALED UP BY 10','SCALED DOWN  BY 10','');
xlabel('BUDGET ITEM NUMBER','FontWeight','bold','FontSize',14);
ylabel('CEP SCALING','FontWeight','bold','FontSize',14);
title('~1 NMi/Hr ''MEDIUM ACCURACY'' INS');
%
% TABULATE RESULTS
%
disp('Error Budget Sensitivities');
disp('Par.                                                  BASE   CEP [NMI/HR]');
disp('No.  ERROR BUDGET PARAMETER NAME AND UNITS            VALUE  /10      x10');
for item=1:7,
    disp([num2str(item),', ',EBname(item,:),'   ',num2str(EB(item)),'   ',num2str(CEPrate1(item)),'   ',num2str(CEPrate2(item))]);
end;
