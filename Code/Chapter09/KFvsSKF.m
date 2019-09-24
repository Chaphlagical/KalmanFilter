% 
disp('M. S. Grewal & A. P. Andrews');
disp('Kalman Filtering Theory and Practice Using MATLAB');
disp('4th Edition, Wiley, 2014.');
disp(' ');
%
% Schmidt-Kalman filtering evaluation using covariance analysis.
% Plots rms estimation uncertainties from solutions of
% covariance equations.
%
disp('Compares Schmidt-Kalman filter vs Kalman filter on problem of');
disp('estimating the state of a damped harmonic resonator excited');
disp('by white noise, and using measurements of resonator displacement');
disp('corrupted by white noise plus exponentially correlated bias.');
disp(' ');
disp('NOTE: This hypothetical example is designed to show large differences');
disp('between Kalman and Schmidt-Kalman filtering.  For most applications of');
disp('Schmidt-Kalman filtering, the penalty for using the suboptimal SKF to');
disp('reduce computation is rather slight.  Nevertheless, it is always best');
disp('to know the increased performace penalty associated with the potentially');
disp('decreased processing costs.');
% 
close all;
clear all;
%
% MODEL PARAMETERS
%
dt     = .002;    % discrete time step (sec)
omega  = 2*pi*30; % damped resonant frequency (30 Hz)
tau(1) = 10/omega;% resonator damping time constant (sec)
tau(2) = tau(1);  % measurement bias correlation time (sec)
p011   = 2^2;     % mean squared resonator displacement (mm^2)
p033   = 1^2;     % steady state mean squared correlated measurement bias (mm^2 equiv)
Rk     = 1^2;     % covariance of white measurement noise
%
%  Dynamic model in discrete time
%
%    Schmidt-Kalman partitioned state transition matrix
%
t1 = 1/tau(1);
t2 = exp(-dt*t1);
t3 = dt*omega;
t4 = cos(t3)*omega*tau(1);
t5 = sin(t3);
t6 = 1/omega;
t7 = t6*t1;
t8 = t2*t5;
t9 = omega^2;
t10 = tau(1)^2;
Phi1(1,1) = t2*(t4+t5)*t7;
Phi1(1,2) = t8*t6;
Phi1(2,1) = -t8*(t9*t10+1)*t6/t10;
Phi1(2,2) = t2*(t4-t5)*t7;
phi33     = exp(-dt/tau(2));
Phi2      = [phi33];
%
%  Kalman model
%
Phik = [Phi1(1,1),Phi1(1,2),0;Phi1(2,1),Phi1(2,2),0;0,0,phi33];
%
% Kalman measurement sensitivity matrix
%
Hk     = [1 0 1]; 
%
% Schmidt-Kalman measurement sensitivity matrices
%
H1     = [1,0];
H2     = [1];
%
%  Stochastic parameters in discrete time
%
%
%    Initialize the Kalman covariance matrix Pk
%    and process noise covariance Qk by solving
%    steady state state covariance equation.
%
Qk = zeros(3);
Pk = zeros(3);
Pk(1,1) = p011;
Pk(3,3) = p033;
t1 = p011;
t2 = dt*omega;
t3 = cos(t2);
t4 = t3*omega;
t5 = tau(1);
t8 = 1/t5;
t9 = dt*t8;
t10 = exp(-2*t9);
t13 = sin(t2);
Pk(1,2) = -t1*(-t4*t5+t5*t3*omega*t10+t13*t10+t13)*t8/t13/(1+t10);
Pk(2,1) = Pk(1,2);
t23 = omega^2;
t25 = exp(2*t9);
t26 = exp(4*t9);
t28 = t3^2;
t30 = t28*t25;
t34 = 1/(-t25+t30-1+t28);
Qk(2,2) = -t23*t1*(-t25+t26-1+exp(6*t9)-4*t28*t26+4*t30)*t10*t34;
t37 = t5^2;
t38 = t23*t37;
t43 = t5*t13;
Pk(2,2) = t1*(-t38*t26-t38*t25+3*t38*t30-t38*t28+2*t43*t4*t25-2*t43*t4-t25+t30-1+t28)/t37*t34;
Qk(3,3) = p033*(1-exp(-2*dt/tau(2)));
%
%  Schmidt-Kalman model process noise covariances
%
Q1      = [0,0;0,Qk(2,2)];
Q2      = [Qk(3,3)];
%
%  Initialize the Schmidt-Kalman covariance matrix
%
P11      = zeros(2);
P11(1,1) = p011;
P11(2,2) = Pk(2,2);
P12      = zeros(2,1);
P22      = [p033];
l        = 0;
%
% initialize plotting arrays
%
sdk1  = 0;
sdk2  = 0;
sdk3  = 0;
sdsk1 = 0;
sdsk2 = 0;
sdsk3 = 0;
t     = 0;
%
  for k=1:50,
  l       = l + 1;
  t(l)    = dt*(k-1); % save a priori values at current time
  sdk1(l) = sqrt(Pk(1,1));
  sdk2(l) = sqrt(Pk(2,2));
  sdk3(l) = sqrt(Pk(3,3));
  sdsk1(l)= sqrt(P11(1,1));
  sdsk2(l)= sqrt(P11(2,2));
  sdsk3(l)= sqrt(P22(1,1));
%
% Kalman observational update
%
  Pk      = Pk - (Pk*Hk'/(Hk*Pk*Hk'+Rk))*Hk*Pk;
  Pk      = .5*(Pk+Pk'); % symmetrize
  l       = l + 1;
  t(l)    = dt*(k-1); % save a posterior values at same time
  sdk1(l) = sqrt(Pk(1,1));
  sdk2(l) = sqrt(Pk(2,2));
  sdk3(l) = sqrt(Pk(3,3));
%
% Schmidt-Kalman gain
%
  C0 = P11*H1'+P12*H2';
  D0 = P12'*H1'+P22*H2';
  E0 = H1*C0+H2*D0;
  Ksk = C0/E0;
%
% Schmidt-Kalman observational update
%
  A0    = eye(2) - Ksk*H1;
  B0    = Ksk*H2;
  P12o  = A0*P12 - B0*P22;
  P11o  = (A0*P11 - B0*P12')*A0' - P12*B0' + Ksk*Rk*Ksk';
  P11   = .5*(P11o+P11o');
  P12   = P12o;
  sdsk1(l)= sqrt(P11(1,1));
  sdsk2(l)= sqrt(P11(2,2));
  sdsk3(l)= sqrt(P22(1,1));
%
% Kalman temporal update
%
  Pk      = Phik*Pk*Phik'+Qk;
  Pk      = .5*(Pk+Pk');
%
% Schmidt-Kalman temporal update
%
  P11   = Phi1*P11*Phi1' + Q1;
  P12   = Phi1*P12*Phi2';
  P22   = Phi2*P22*Phi2' + Q2;
  P11   = .5*(P11+P11');
  P22   = .5*(P22+P22');
  end;
subplot(3,1,1),plot(t,sdk1,'r-',t,sdsk1,'b-');ylabel('\sigma_{position}');
title('Resonator tracking using measurements with exponentially correlated bias');
subplot(3,1,2),plot(t,sdk2,'r-',t,sdsk2,'b-');ylabel('\sigma_{velocity}');
subplot(3,1,3),plot(t,sdk3,'r-',t,sdsk3,'b-');
hold on;
plot(.9*max(t),1.2*max(sdsk3),'w.');
plot(.9*max(t),.1*max(sdsk3),'w.');
hold off;
text(.9*max(t),1.05*max(sdsk3),'SKF','Color','b','VerticalAlignment','bottom','HorizontalAlignment','right');
text(.7*max(t),0.95*min(sdk3),'KF','Color','r','VerticalAlignment','top','HorizontalAlignment','right');
ylabel('\sigma_{meas. bias}');
xlabel('Time');


