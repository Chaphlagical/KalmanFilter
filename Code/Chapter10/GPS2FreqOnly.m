%%   
%%  Grewal & Andrews,
%%  Kalman Filtering: Theory and Practice Usin MATLAB, 4th Edition,
%%  Wiley, 2014.
%%  
%
% PERFORMANCE ANALYSIS OF DUAL-FREQUENCY GPS RECEIVER NAVIGATING ON A
% 100-KM FIGURE-8 TEST TRACK, WITH SIGNAL OUTAGES 1 MINUTE OUT OF EVERY HOUR
%
% STATE VECTOR:
%   9 HOST VEHICLE MODEL DYNAMIC VARIABLES
%   2 RECEIVER CLOCK ERROR VERIABLES
%
close all;
clear all;
%
% Earth model
%
REarth        = 6371009;    % mean radius of earth [m]
%
% Test site location
%
Lat       = 39.7931*pi/180;         % 
Lon       = -86.2389*pi/180;        % Approximate location of
Alt       = 221;                    % Indianapolis Motor Speedway
LatRad    = Lat*pi/180;             %
LonRad    = Lon*pi/189;             %
cLat      = cos(LatRad);
%
% 15-degree horizon mask limit for satellite selection
%
s15 = sin(15*pi/180); % sine of mask angle
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
%   NINE HOST VEHICLE DYNAMIC MODEL STATE VARIABLES
%
%   1. xENU(1)  Easting error [m]
%   2. xENU(2)  Northing error [m]
%   3. xENU(3)  Altitude error [m]
%   4. vENU(1)  East velocity error [m/s]
%   5. vENU(2)  North velocity error [m/s]
%   6. vENU(3)  Altitude rate error [m/s]
%   7. aENU(1)  East acceleration error [m]
%   8. aENU(2)  North acceleration error [m]
%   9. aENU(3)  Upward acceleration error [m]
%
%  HOST VEHICLE DYNAMIC MODEL PARAMETERS ARE "TUNED" USING
%  STATISTICS FROM Big8TrackStats.m (SIMULATOR FOR 100KM FIGURE-8 TRACK)
%  
%  RMS E-W Position Excursion   = 4566.5902 meter
%  RMS N-S Position Excursion   = 13699.7706 meter
%  RMS Vert. Position Excursion = 6.1237 meter
%  RMS E-W Velocity             = 15.9426 m/s
%  RMS N-S Velocity             = 23.9139 m/s
%  RMS Vert. Velocity           = 0.0061707 m/s
%  RMS E-W Acceleration         = 0.055643 m/s/s
%  RMS N-S Acceleration         = 0.041732 m/s/s
%  RMS Vert. Acceleration       = 1.0771e-005 m/s/s
%  RMS Delta Velocity East      = 0.055642 m/s at Delta t = 1 sec.
%  RMS Delta Velocity North     = 0.041732 m/s at Delta t = 1 sec.
%  RMS Delta Velocity Up        = 1.077e-005 m/s at Delta t = 1 sec.
%  E. Position Correlation Time = 354.155 sec
%  N. Position Correlation Time = 735.6834 sec
%  Vertical Position Corr. Time = 1661.1583 sec
%  E. Velocity Correlation Time = 330.2384 sec
%  N. Velocity Correlation Time = 635.1474 sec
%  Vertical Velocity Corr. Time = 735.6834 sec
%  E. Acceler. Correlation Time = 354.155 sec
%  N. Acceler. Correlation Time = 735.6834 sec
%  Vertical Acceler. Corr. Time = 635.1474 sec
%
% "TUNING" USES THESE STATISTICS IN THE DAMP3 MODELS FOR EACH AXIS
% (East, North ,and Up), then re-sorts the resulting parameter matrices 
% into the order of the state variables shown above.
%
% The tuning is accomplished within the callable function
%
% [Phi,Q,P] = DAMP3ParamsD(SigmaAcc,SigmaVel,SigmaPos,TauAcc,DeltaT)
%
% by plugging in the appropriate statistical values from above.
%
[PhiE,QE,PE] = DAMP3ParamsD(0.055643,15.9426,4566.5902,354.155,Dt);
[PhiN,QN,PN] = DAMP3ParamsD(0.041732,23.9139,13699.7706,735.6834,Dt);
[PhiU,QU,PU] = DAMP3ParamsD(1.0771e-5,0.0061707,6.1237,635.1474,Dt);
%
% However, re-arranging the returned matrix parameters into state-vector
% order requires they be transformed from three independent 3-state filters
% for East, North, and Up to state-vector coordinate.
%
% The 9x9 transformation matrix from Pos-Vel-Acc coordinates to ENU
% coordinates is a symmetric orthogonal matrix such that inv(T) == T' == T:
%
T   = zeros(9);
T(1,1) = 1;
T(2,4) = 1;
T(3,7) = 1;
T(4,2) = 1;
T(5,5) = 1;
T(6,8) = 1;
T(7,3) = 1;
T(8,6) = 1;
T(9,9) = 1;
%
% The resulting transformed model matrices are then
%
PhiVeh  = T*[PhiE,zeros(3,6);zeros(3),PhiN,zeros(3);zeros(3,6),PhiU]*T';
PVeh    = T*[PE,zeros(3,6);zeros(3),PN,zeros(3);zeros(3,6),PU]*T';
QVeh    = T*[QE,zeros(3,6);zeros(3),QN,zeros(3);zeros(3,6),QU]*T';
%
%
%
Phi     = [PhiVeh,zeros(9,2);zeros(2,9),Phiclock];
P       = [PVeh,zeros(9,2);zeros(2,9),P0clock];
P       = .5*(P+P');
Q       = [QVeh,zeros(9,2);zeros(2,9),Qclock*Dt];
Q       = .5*(Q+Q');
%
% Initialize plotting variable arrays
%
CEPGNSS     = [];
PosGNSS     = [];
VelGNSS     = [];
AccGNSS     = [];
NoSatsUsed  = []; % Number of satellites available for navigation
time        = []; % Simulation time [s] 
RMSVert     = [];
RMSclock    = []; % RMS clock timing uncertainty [m]
RMSdrift    = []; % RMS clock frequency uncertainty [m/s]
RMSHorGPS   = [];
RMStrue     = 0;
RMSest      = 0;
Nrms        = 0;
EposRMS     = []; % RMS East position uncertainty
NposRMS     = []; % RMS North position uncertainty
EvelRMS     = []; % RMS East velocity uncertainty
NvelRMS     = []; % RMS North velocity uncertainty
%
% Set test length
%
hours       = 4;
minutes     = hours*60;
seconds     = minutes*60;
%
%
for t=0:Dt:seconds,
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
    H           = [HGNSS,zeros(NumSats,6),ones(NumSats,1),zeros(NumSats,1)];
    beforetunnel = (mod(t,3600) < (3600 - 120));
    aftertunnel  = (mod(t,3600) > (3600 - 60));
    if beforetunnel|aftertunnel,
        PHT         = P*H';
        Denom       = H*PHT + R;
        K           = PHT/Denom;
        P           = P - K*PHT';
        P           = .5*(P + P');
    else
        NoSatsAvail = 0;
    end;
    time        = [time,t];
    NoSatsUsed  = [NoSatsUsed,NoSatsAvail];
    EposRMS     = [EposRMS,sqrt(P(1,1))]; % RMS East position uncertainty
    NposRMS     = [NposRMS,sqrt(P(2,2))]; % RMS North position uncertainty
    EvelRMS     = [EvelRMS,sqrt(P(4,4))]; % RMS East velocity uncertainty
    NvelRMS     = [NvelRMS,sqrt(P(5,5))]; % RMS North velocity uncertainty
    RMSHorGPS   = [RMSHorGPS,sqrt(P(1,1)+P(2,2))];
    CEPGNSS     = [CEPGNSS,1.2*sqrt(P(1,1)+P(2,2))/1852];
    PosGNSS     = [PosGNSS,sqrt(P(1,1)+P(2,2)+P(3,3))];
    VelGNSS     = [VelGNSS,sqrt(P(4,4)+P(5,5)+P(6,6))];
    AccGNSS     = [AccGNSS,sqrt(P(7,7)+P(8,8)+P(9,9))];
    RMSVert     = [RMSVert,sqrt(P(3,3))]; 
    RMSclock    = [RMSclock,sqrt(P(10,10))];
    RMSdrift    = [RMSdrift,sqrt(P(11,11))];
    P           = Phi*P*Phi' + Q;
    time        = [time,t];
    NoSatsUsed  = [NoSatsUsed,NoSatsAvail];
    EposRMS     = [EposRMS,sqrt(P(1,1))]; % RMS East position uncertainty
    NposRMS     = [NposRMS,sqrt(P(2,2))]; % RMS North position uncertainty
    EvelRMS     = [EvelRMS,sqrt(P(4,4))]; % RMS East velocity uncertainty
    NvelRMS     = [NvelRMS,sqrt(P(5,5))]; % RMS North velocity uncertainty
    RMSHorGPS   = [RMSHorGPS,sqrt(P(1,1)+P(2,2))];
    CEPGNSS     = [CEPGNSS,1.2*sqrt(P(1,1)+P(2,2))/1852];
    PosGNSS     = [PosGNSS,sqrt(P(1,1)+P(2,2)+P(3,3))];
    VelGNSS     = [VelGNSS,sqrt(P(4,4)+P(5,5)+P(6,6))];
    AccGNSS     = [AccGNSS,sqrt(P(7,7)+P(8,8)+P(9,9))];
    RMSVert     = [RMSVert,sqrt(P(3,3))]; 
    RMSclock    = [RMSclock,sqrt(P(10,10))];
    RMSdrift    = [RMSdrift,sqrt(P(11,11))];
end;
save 'RMSHorGPS2' 'RMSHorGPS'
endoftime = length(time); % (100:endoftime)
figure;
plot(time(100:endoftime)/3600,RMSHorGPS(100:endoftime),'k-','LineWidth',2);
set(gca,'FontSize',14);
title('DUAL-FREQUENCY GPS RECEIVER','FontWeight','bold','FontSize',14);
xlabel('TIME [hr]');
ylabel('RMS HORIZONTAL UNCERT. [m]','FontWeight','bold','FontSize',14);
figure;
subplot(2,2,1),
plot(time(100:endoftime)/3600,CEPGNSS(100:endoftime)*1852,'k-',time(100:endoftime)/3600,RMSVert(100:endoftime),'k--','LineWidth',2);
set(gca,'FontSize',14);
title('GNSS ON 100-KM FIGURE-8 TRACK','FontWeight','bold','FontSize',12);
xlabel('TIME [HR]','FontWeight','bold','FontSize',12);
ylabel('RMS POS. [m]','FontWeight','bold','FontSize',12);
%legend('Hor.','Vert.');
subplot(2,2,2),
plot(time(100:endoftime)/3600,VelGNSS(100:endoftime),'k-','LineWidth',2);
set(gca,'FontSize',14);
xlabel('TIME [HR]','FontWeight','bold','FontSize',12);
ylabel('RMS VEL. [m/s]','FontWeight','bold','FontSize',12);
subplot(2,2,3),
plot(time(100:endoftime)/3600,AccGNSS(100:endoftime),'k-','LineWidth',2);
set(gca,'FontSize',14);
xlabel('TIME [HR]','FontWeight','bold','FontSize',12);
ylabel('RMS ACC [m/s/s]','FontWeight','bold','FontSize',12);
subplot(2,2,4),
plot(time/3600,NoSatsUsed,'k-','LineWidth',2);
set(gca,'FontSize',14);
xlabel('TIME [HR]','FontWeight','bold','FontSize',12);
ylabel('NO. SATs USED','FontWeight','bold','FontSize',12);
figure;
subplot(2,1,1),
plot(time(100:endoftime)/3600,RMSclock(100:endoftime),'k-','LineWidth',2);
set(gca,'FontSize',14);
title('GNSS ON 100-KM FIGURE-8 TRACK','FontWeight','bold','FontSize',12);
xlabel('TIME [HR]','FontWeight','bold','FontSize',12);
ylabel('RMS CLOCK ERR. [m]','FontWeight','bold','FontSize',12);
subplot(2,1,2),
plot(time(100:endoftime)/3600,RMSdrift(100:endoftime),'k-','LineWidth',2);
set(gca,'FontSize',14);
xlabel('TIME [HR]','FontWeight','bold','FontSize',12);
ylabel('CLOCK DRIFT [m/s]','FontWeight','bold','FontSize',12);
figure;
%
EoRMS = length(EposRMS);  % (100:EoRMS)
subplot(2,1,1),
plot(EposRMS(100:EoRMS),'k-','LineWidth',2);
set(gca,'FontSize',14);
xlabel('SAMPLE NUMBER');
ylabel('RMS EAST POS. ERR. [m]','FontWeight','bold','FontSize',12);
title('DUAL-FREQUENCY GPS RECEIVER','FontWeight','bold','FontSize',12);
subplot(2,1,2),
plot(NposRMS(100:EoRMS),'k-','LineWidth',2);
set(gca,'FontSize',14);
xlabel('SAMPLE NUMBER');
ylabel('RMS NORTH POS. ERR. [m]','FontWeight','bold','FontSize',12);
figure;
subplot(2,1,1),
plot(EvelRMS(100:EoRMS),'k-','LineWidth',2);
set(gca,'FontSize',14);
xlabel('SAMPLE NUMBER');
ylabel('RMS E. VEL. ERR. [m/s]','FontWeight','bold','FontSize',12);
title('GNSS ON 100-KM FIGURE-8 TRACK','FontWeight','bold','FontSize',12);
subplot(2,1,2),
plot(NvelRMS(100:EoRMS),'k-','LineWidth',2);
set(gca,'FontSize',14);
ylabel('RMS N. VEL. ERR [m/s]','FontWeight','bold','FontSize',12);
xlabel('SAMPLE NUMBER');
