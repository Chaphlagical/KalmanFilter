%
% This demonstration uses the problem in Example 4.3
% to demonstrate the effects of
% approximations commonly used in converting
% continuous time models to discrete time models.
%
close all;
clear all;
disp('Grewal & Andrews,');
disp('Kalman Filtering: Theory and Practice Using MATLAB, 4th Edition');
disp('Wiley, 2014.');
disp(' ');
disp('Demonstration of errors due to first-order approximations in');
disp('converting a continuous-time model with dynamic coefficient matrix F');
disp('to a discrete-time model with a state transition matrix Phi.');
disp(' ');
%
% Approximations include the following three:
%
% 1   Approximating the state transition matrix (Phi)
%     by the first two terms in the power series
%     expansion of Phi = exp(T*F), which is
%
%        Phia = I + T*F,
%
%     where
%
%       Phia is the approximated value of Phi
%       I    is a 2x2 identity matrix
%       T    is the sample interval in seconds
%       F    is the dynamic coefficient matrix
%
% 2   Approximating the forcing step in discrete time
%     as delta = T*f, where f is the forcing function
%     in the differential equation x-dot = F x + f and
%     delta is the additive forcing step:
%        x(k+1) = Phi * x(k) + delta
%     in the discrete time model.
%
% 3   Approximating the additive process noise covariance
%     as Qd = T* Qc, where Qd is the process noise
%     covariance the discrete time Riccati equation 
%        P(k+1) = Phi*P(k)*Phi^T + Qd
%     and Qc is the covariance in the continuous equation
%        P-dot = F*P + P*F^T + Qc
%     The exact solution is defined by an integral
%                      T
%        Qd = Phi*[ int  e^(-sF)*Qc*e^(-sF^T) ds ]*Phi^T
%                      0
%
% The parameters of the problem are:
%
I     = eye(2);% (identity matrix)
T     = .02;   % intersample time interval (s)
               % (not specified in Example 4.3)
qc    = 4.47;  % continuous process noise covariance (ft^2/s^3)
R     = .01;   % measurement noise variance (s*ft^2)
zeta  = 0.2;   % damping factor (unitless)
omega = 5;     % undamped resonant frequency (rad/s)
f     = [0;12];% non-homogeneous forcing (constant)
H     = [1,0]; % measurement sensitivity matrix
%
% The dynamic coefficient matrix is a function of these:
%
F      = [0,1;-omega^2,-2*zeta*omega];
%
% Approximated discrete time model parameters
%
Phia   = I + T*F; % approximated Phi
deltaa = T*f;     % approximated forcing step
Qda    = [0,0;0,qc*T];
%
% Exact values of the model parameters:
%  (using intermediate parameters for efficiency)
%
lambda = exp(-T*omega*T);
psi    = 1 - zeta^2;
xi     = sqrt(psi);
theta  = xi*omega*T;
c      = cos(theta);
s      = sin(theta);
Phi    = zeros(2);
%
% The exact solution for exp(T*F) is
%
Phi(1,1) = lambda*c + zeta*s/xi;
Phi(1,2) = lambda*s/(omega*xi);
Phi(2,1) = -omega*lambda*s/xi;
Phi(2,2) = lambda*c - zeta*s/xi;
%
% and the exact forcing term is (check this)
%
delta = f(2)*[(1 - lambda*(c-zeta*s*xi/psi))/omega^2;
              lambda*s/(omega*xi)];
%
% discrete time process noise covariance
%
l1  = zeta^2;
l4  = exp(-2*omega*zeta*T);
l6  = sqrt(1-l1);
l8  = l6*omega*T;
l9  = cos(2*l8);
l11 = l4*l9*l1;
l15 = l4*sin(2*l8)*l6*zeta;
l19 = 1/(l1-1);
l20 = omega^2;
l24 = 1/zeta;
l32 = l4*(l9-1)*l19/l20;
Qd  = (qc/4)*[(l1-l11+l15-1+l4)*l19*l24,l32;
             l32,(l1-l11-l15-1+l4)*l19*l24/omega];
disp(' ');
disp('Exact and approximated values of state transition matrix (Phi)');
disp([Phi,Phia]);
disp('Exact and approximated values of forcing step (delta)');
disp([delta,deltaa]);
disp('Exact and approximated values of process noise cov. (Qd)');
disp([Qd,Qda]);
disp('(See the script for details.)');
disp(' ');
disp('(Allow a moment for simulation.)');
%
% Initial values:
%
x    = [0;0];       % "true" state value
xa   = [0;0];       % approximated value
P    = [2,0;0,2];   % covariance of estimation uncertainty
Pa   = [2,0;0,2];   % ditto, approximated
dx   = [sqrt(P(1,1))*randn(1);sqrt(P(2,2))*randn(1)];
xh   = x + dx;  % estimated state
xha  = xa + dx; % approximated state
%
% The following are for pseudonoise scaling
%
RMSm = sqrt(R);     % RMS measurement noise
RMSn = zeros(2);    % Process noise weighting matrix...
RMSn(2,2) = sqrt(Qd(2,2)); % is upper triangular...
RMSn(1,2) = Qd(1,2)/RMSn(2,2); % cholesky factor...
RMSn(1,1) = sqrt(Qd(1,1)-RMSn(1,2)^2); % of Qd.
%
% Initialize arrays for simulated values
%
m=0;tm=0;t=0;x1=0;x2=0;xa1=0;xa2=0;K1=0;K2=0;
z=0;xh1=0;xh2=0;u1=0;l1=0;u2=0;l2=0;sn1=0;sn2=0;
%
% Simulate and save values (identified by comments)
%
   for k=1:101,
   t(k)   = (k-1)*T; % time starting with t=0
   x1(k)  = x(1);    % true state variables
   x2(k)  = x(2);
   xa1(k) = xa(1);   % approximated state variables
   xa2(k) = xa(2);
   noise  = RMSn*[randn(1);randn(1)]; % process noise
   x      = Phi*x + delta + noise;
   xa     = Phia*xa + deltaa + noise;
   zn     = x(1) + RMSm*randn(1); % Measurement equals true
   z(k)   = zn;                    % position plus noise.
%
% The forcing term has no effect on P.  Why?
%
   xh    = Phi*xh + delta;        % a priori estimates
   xha   = Phia*xha + deltaa;
%
% Use double indexing for priori/posteriori values
%
   m      = m+1;
   tm(m)  = t(k);  % time of a priori value
   xh1(m) = xh(1); % a priori estimates
   xh2(m) = xh(2);
   xha1(m)= xha(1); % (approx)
   xha2(m)= xha(2);
   P      = Phi*P*Phi' + Qd;
   K      = P*H'/(H*P*H'+R);
   K1(k)  = K(1); % Kalman gains
   K2(k)  = K(2);
   Pa     = Phia*Pa*Phia' + Qda;   % Approx. temporal update
   Ka     = Pa*H'/(H*Pa*H'+R);
   Ka1(k) = Ka(1); % approximated Kalman gains
   Ka2(k) = Ka(2);
   s1(m)  = sqrt(P(1,1)); % RMS uncertainties of estimates
   s2(m)  = sqrt(P(2,2)); % (exact)
   sa1(m) = sqrt(Pa(1,1)); % RMS uncertainties of estimates
   sa2(m) = sqrt(Pa(2,2)); % (approx)
%
% posteriori values
%
   m     = m+1;
   tm(m) = t(k);  % Time of a posteriori values.
   xh     = xh + K*(zn - H*xh); % a posteriori estimate
   xha    = xha+ Ka*(zn - H*xha); % a posteriori estimate
   xh1(m) = xh(1); % a posteriori estimates
   xh2(m) = xh(2);
   xha1(m)= xha(1); % (approx)
   xha2(m)= xha(2);
   P     = P - K*H*P; 
   P     = .5*(P+P'); % preserve symmetry of P
   Pa    = Pa - Ka*H*Pa;
   Pa   = .5*(Pa+Pa');
   s1(m)  = sqrt(P(1,1)); % RMS uncertainties of estimates
   s2(m)  = sqrt(P(2,2)); % (exact)
   sa1(m) = sqrt(Pa(1,1)); % RMS uncertainties of estimates
   sa2(m) = sqrt(Pa(2,2)); % (approx)
   xh1(m)= xh(1); % a posteriori estimates
   xh2(m)= xh(2);
   end;
disp('Done.');
disp(' ');
disp('Behold the plots in the accompanying figure.');
disp(' ');
disp('These plots show the "true" state dynamics');
disp('and values using approximated parameters.');
disp(' ');
disp('Is the error due to approximation discernable?');
disp(' ');
disp('Press <ENTER> key to continue');
figure;
subplot(2,1,1),plot(t,x1,'b-',t,xa1,'k-');xlabel('Time (sec)');ylabel('Position');title('F to \Phi Approximation: true (blue) and approx. (blk) state dynamics');
subplot(2,1,2),plot(t,x2,'b-',t,xa2,'k-');xlabel('Time (sec)');ylabel('Velocity');
pause
%
% The second set of plots are the rms estimation uncertainties
% and correlation coefficient as a function of time for the
% the two cases (exact and approximated parameters).
%
disp(' ');
disp('Displayed plots show RMS state uncertainties');
disp('for the exact and approximated models.');
disp(' ');
disp('Press <ENTER> key to continue');
figure;
subplot(2,1,1),plot(tm,s1,'b',tm,sa1,'k-');xlabel('Time (sec)');ylabel('RMS Pos.');title('F to \Phi Approximation: (uncertainties) true (blue) and approx. (blk)');
subplot(2,1,2),plot(tm,s2,'b',tm,sa2,'k-');xlabel('Time (sec)');ylabel('RMS Vel.');
pause
disp(' ');
disp('Displayed state estimation plots are:');
disp(' 1. True, estimated and measured position (top).');
disp('    (including estimates with exact and approx. models)');
disp(' 2. True and estimated velocity (bottom).');
disp('    (including estimates with exact and approx. models)');
disp(' ');
disp('NOTE: Different pseudorandom sequences will make these');
disp('      plots different on repeated runs.');
disp(' ');
figure;
subplot(2,1,1),plot(t,x1,'b-',t,z,'r-',tm,xh1,'g-',tm,xha1,'k-');xlabel('Time (sec)');ylabel('Position');title('F to \Phi Approximation: true (blu), estim. (grn) est/approx (blk) & meas. (red) states');
subplot(2,1,2),plot(t,x2,'b-',tm,xh2,'g-',tm,xha2,'k-');xlabel('Time (sec)');ylabel('Velocity');
