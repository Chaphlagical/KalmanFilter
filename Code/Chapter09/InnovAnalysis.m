%
% Grewal & Andrews,
% Kalman Filtering: Theory and Practice Using MATLAB, 4th Edition,
% Wiley, 2014.
%
% Example of innovations analysis with model parameters
% scaled up and down by a factor of two, to demonstrate
% its effect on weighted squared innovations.
%
% The model is a scalar exponentially correlated process.
%
close all;
clear all;
%
DeltaT = 1;         % Time-step
tau    = 10*DeltaT; % Correlation time.
%
% "True" model
%
Phi0   = exp(-DeltaT/tau);
P0     = 1;  % steady-state signal
Q0     = P0*(1-Phi0^2);
R0     = P0/10; % 10dB SNR
H0     = 1;
x0     = 0;
xtrue  = sqrt(P0)*randn;
%
% Models with Phi scaled up and down by a factor of two
%
x1     = x0;
x2     = x0;
P1     = P0;
P2     = P0;
Phi1   = Phi0*2;
Phi2   = Phi0/2;
R1     = R0;
R2     = R0;
Q1     = Q0;
Q2     = Q0;
H1     = H0;
H2     = H0;
%
% Models with H scaled up and down by a factor of two
%
x3     = x0;
x4     = x0;
P3     = P0;
P4     = P0;
Phi3   = Phi0;
Phi4   = Phi0;
R3     = R0;
R4     = R0;
Q3     = Q0;
Q4     = Q0;
H3     = H0*2;
H4     = H0/2;
%
% Models with Q scaled up and down by a factor of two
%
x5     = x0;
x6     = x0;
P5     = P0;
P6     = P0;
Phi5   = Phi0;
Phi6   = Phi0;
R5     = R0;
R6     = R0;
Q5     = Q0*2;
Q6     = Q0/2;
H5     = H0;
H6     = H0;
%
% Models with R scaled up and down by a factor of two
%
x7     = x0;
x8     = x0;
P7     = P0;
P8     = P0;
Phi7   = Phi0;
Phi8   = Phi0;
R7     = R0*2;
R8     = R0/2;
Q7     = Q0;
Q8     = Q0;
H7     = H0;
H8     = H0;
Prior0 = 0;
Post0  = 0;
Prior1 = 0;
Post1  = 0;
Prior2 = 0;
Post2  = 0;
Prior3 = 0;
Post3  = 0;
Prior4 = 0;
Post4  = 0;
Prior5 = 0;
Post5  = 0;
Prior6 = 0;
Post6  = 0;
Prior7 = 0;
Post7  = 0;
Prior8 = 0;
Post8  = 0;
m      = 0;
Nk     = 1001;
for k=1:Nk,
    m       = m + 1;
    Tm(m)   = k-1;
    Tk(k)   = k-1;
    XT(m)   = xtrue;
    X0(m)   = x0;
    Prior0  = Prior0 + (x0 - xtrue)^2;
    z       = H0*xtrue+sqrt(R0)*randn;
    dz0     = z - H0*x0;
    %
    % Measurement update
    %
    PH0     = P0*H0;
    Yzz0    = 1/(H0*PH0+R0);
    Chi0(k) = dz0^2*Yzz0;
    K0      = PH0*Yzz0;
    x0      = x0 + K0*dz0;
    Post0   = Post0 + (x0 - xtrue)^2;
    P0      = P0 - K0*PH0;
    m       = m + 1;
    Tm(m)   = k-1;
    XT(m)   = xtrue;
    X0(m)   = x0;
    %
    % Time update
    %
    xtrue   = Phi0*xtrue + sqrt(Q0)*randn;
    x0      = Phi0*x0;
    P0      = Phi0^2*P0 + Q0;
    %
    % Repeat for other models
    %
    Prior1  = Prior1 + (x1 - xtrue)^2;
    Prior2  = Prior2 + (x2 - xtrue)^2;
    Prior3  = Prior3 + (x3 - xtrue)^2;
    Prior4  = Prior4 + (x4 - xtrue)^2;
    Prior5  = Prior5 + (x5 - xtrue)^2;
    Prior6  = Prior6 + (x6 - xtrue)^2;
    Prior7  = Prior7 + (x7 - xtrue)^2;
    Prior8  = Prior8 + (x8 - xtrue)^2;
    PH1     = P1*H1;
    PH2     = P2*H2;
    PH3     = P3*H3;
    PH4     = P4*H4;
    PH5     = P5*H5;
    PH6     = P6*H6;
    PH7     = P7*H7;
    PH8     = P8*H8;
    Yzz1    = 1/(H1*PH1+R1);
    Yzz2    = 1/(H2*PH2+R2);
    Yzz3    = 1/(H3*PH3+R3);
    Yzz4    = 1/(H4*PH4+R4);
    Yzz5    = 1/(H5*PH5+R5);
    Yzz6    = 1/(H6*PH6+R6);
    Yzz7    = 1/(H7*PH7+R7);
    Yzz8    = 1/(H8*PH8+R8);
    dz1     = z - H1*x1;
    dz2     = z - H2*x2;
    dz3     = z - H3*x3;
    dz4     = z - H4*x4;
    dz5     = z - H5*x5;
    dz6     = z - H6*x6;
    dz7     = z - H7*x7;
    dz8     = z - H8*x8;
    Chi1(k) = dz1^2*Yzz1;
    Chi2(k) = dz2^2*Yzz2;
    Chi3(k) = dz3^2*Yzz3;
    Chi4(k) = dz4^2*Yzz4;
    Chi5(k) = dz5^2*Yzz5;
    Chi6(k) = dz6^2*Yzz6;
    Chi7(k) = dz7^2*Yzz7;
    Chi8(k) = dz8^2*Yzz8;
    K1      = PH1*Yzz1;
    K2      = PH2*Yzz2;
    K3      = PH3*Yzz3;
    K4      = PH4*Yzz4;
    K5      = PH5*Yzz5;
    K6      = PH6*Yzz6;
    K7      = PH7*Yzz7;
    K8      = PH8*Yzz8;
    x1      = x1 + K1*dz1;
    x2      = x2 + K2*dz2;
    x3      = x3 + K3*dz3;
    x4      = x4 + K4*dz4;
    x5      = x5 + K5*dz5;
    x6      = x6 + K6*dz6;
    x7      = x7 + K7*dz7;
    x8      = x8 + K8*dz8;
    Post1   = Post1 + (x1 - xtrue)^2;
    Post2   = Post2 + (x2 - xtrue)^2;
    Post3   = Post3 + (x3 - xtrue)^2;
    Post4   = Post4 + (x4 - xtrue)^2;
    Post5   = Post5 + (x5 - xtrue)^2;
    Post6   = Post6 + (x6 - xtrue)^2;
    Post7   = Post7 + (x7 - xtrue)^2;
    Post8   = Post8 + (x8 - xtrue)^2;
    P1      = P1 - K1*PH1;
    P2      = P2 - K2*PH2;
    P3      = P3 - K3*PH3;
    P4      = P4 - K4*PH4;
    P5      = P5 - K5*PH5;
    P6      = P6 - K6*PH6;
    P7      = P7 - K7*PH7;
    P8      = P8 - K8*PH8;
    %
    % Time update
    %
    x1      = Phi1*x1;
    x2      = Phi2*x2;
    x3      = Phi3*x3;
    x4      = Phi4*x4;
    x5      = Phi5*x5;
    x6      = Phi6*x6;
    x7      = Phi7*x7;
    x8      = Phi8*x8;
    P1      = Phi1^2*P1 + Q1;
    P2      = Phi2^2*P2 + Q2;
    P3      = Phi3^2*P3 + Q3;
    P4      = Phi4^2*P4 + Q4;
    P5      = Phi5^2*P5 + Q5;
    P6      = Phi6^2*P6 + Q6;
    P7      = Phi7^2*P7 + Q7;
    P8      = Phi8^2*P8 + Q8;
end;
%
RMS0    = mean(Chi0);
RMS1    = mean(Chi1);
RMS2    = mean(Chi2);
RMS3    = mean(Chi3);
RMS4    = mean(Chi4);
RMS5    = mean(Chi5);
RMS6    = mean(Chi6);
RMS7    = mean(Chi7);
RMS8    = mean(Chi8);
Prior0  = sqrt(Prior0/Nk);
Post0   = sqrt(Post0/Nk);
Prior1  = sqrt(Prior1/Nk);
Post1   = sqrt(Post1/Nk);
Prior2  = sqrt(Prior2/Nk);
Post2   = sqrt(Post2/Nk);
Prior3  = sqrt(Prior3/Nk);
Post3   = sqrt(Post3/Nk);
Prior4  = sqrt(Prior4/Nk);
Post4   = sqrt(Post4/Nk);
Prior5  = sqrt(Prior5/Nk);
Post5   = sqrt(Post5/Nk);
Prior6  = sqrt(Prior6/Nk);
Post6   = sqrt(Post6/Nk);
Prior7  = sqrt(Prior7/Nk);
Post7   = sqrt(Post7/Nk);
Prior8  = sqrt(Prior8/Nk);
Post8   = sqrt(Post8/Nk);
disp(['$\Phi$ & $H$ & $Q$ & $R$ & ',num2str(Prior0),' & ',num2str(Post0),' & ',num2str(RMS0)]);
disp(['$2\Phi$ & $H$ & $Q$ & $R$ & ',num2str(Prior1),' & ',num2str(Post1),' & ',num2str(RMS1)]);
disp(['$\Phi/2$ & $H$ & $Q$ & $R$ & ',num2str(Prior2),' & ',num2str(Post2),' & ',num2str(RMS2)]);
disp(['$\Phi$ & $2H$ & $Q$ & $R$ & ',num2str(Prior3),' & ',num2str(Post3),' & ',num2str(RMS3)]);
disp(['$\Phi$ & $H/2$ & $Q$ & $R$ & ',num2str(Prior4),' & ',num2str(Post4),' & ',num2str(RMS4)]);
disp(['$\Phi$ & $H$ & $2Q$ & $R$ & ',num2str(Prior5),' & ',num2str(Post5),' & ',num2str(RMS5)]);
disp(['$\Phi$ & $H$ & $Q/2$ & $R$ & ',num2str(Prior6),' & ',num2str(Post6),' & ',num2str(RMS6)]);
disp(['$\Phi$ & $H$ & $Q$ & $2R$ & ',num2str(Prior7),' & ',num2str(Post7),' & ',num2str(RMS7)]);
disp(['$\Phi$ & $H$ & $Q$ & $R/2$ & ',num2str(Prior8),' & ',num2str(Post8),' & ',num2str(RMS8)]);
%
figure;
plot(Tm,XT,'b-',Tm,X0,'r-');
legend('true','est.');
title('Correctly Modeled');
xlabel('Time [sec]');
ylabel('State Value');
%
figure;
subplot(3,1,1);
plot(Tk,Chi0,'k-');
title(['Correct Model: RMS = ',num2str(RMS0)]);
ylabel('{\chi}^2');
subplot(3,1,2);
plot(Tk,Chi1,'k-');
title(['{\Phi}{\times}2: RMS = ',num2str(RMS1)]);
ylabel('{\chi}^2');
subplot(3,1,3);
plot(Tk,Chi2,'k-');
title(['{\Phi}/2: RMS = ',num2str(RMS2)]);
ylabel('{\chi}^2');
xlabel('Time [sec]');
%
figure;
subplot(3,1,1);
plot(Tk,Chi0,'k-');
title(['Correct Model: RMS = ',num2str(RMS0)]);
ylabel('{\chi}^2');
subplot(3,1,2);
plot(Tk,Chi3,'k-');
title(['H{\times}2: RMS = ',num2str(RMS3)]);
ylabel('{\chi}^2');
subplot(3,1,3);
plot(Tk,Chi4,'k-');
title(['H/2: RMS = ',num2str(RMS4)]);
ylabel('{\chi}^2');
xlabel('Time [sec]');
%
figure;
subplot(3,1,1);
plot(Tk,Chi0,'k-');
title(['Correct Model: RMS = ',num2str(RMS0)]);
ylabel('{\chi}^2');
subplot(3,1,2);
plot(Tk,Chi5,'k-');
title(['Q{\times}2: RMS = ',num2str(RMS5)]);
ylabel('{\chi}^2');
subplot(3,1,3);
plot(Tk,Chi6,'k-');
title(['Q/2: RMS = ',num2str(RMS6)]);
ylabel('{\chi}^2');
xlabel('Time [sec]');
%
figure;
subplot(3,1,1);
plot(Tk,Chi0,'k-');
title(['Correct Model: RMS = ',num2str(RMS0)]);
ylabel('{\chi}^2');
subplot(3,1,2);
plot(Tk,Chi7,'k-');
title(['R{\times}2: RMS = ',num2str(RMS7)]);
ylabel('{\chi}^2');
subplot(3,1,3);
plot(Tk,Chi8,'k-');
title(['R/2: RMS = ',num2str(RMS8)]);
ylabel('{\chi}^2');
xlabel('Time [sec]');
%
