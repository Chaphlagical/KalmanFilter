disp('Grewal & Andrews');
disp('Kalman Filtering: Theory and Practice Using MATLAB, 4th Edition,');
disp('Wiley, 2014');
disp(' ');
disp('Example of 3-pass fixed-interval smoothing, with a unit-variance');
disp('exponentially correlated dynamic model.');
disp(' ');
%
close all;
clear all;
%
% Model parameters
%
Pinf = 1;                       % variance of process (w/ no measurements)
tau  = 10;                       % correlation time-constant
Dt   = 1;                       % discrete time step
H    = 1;                       % measurement sensitivity
F    = -1/tau;                  % dynamic coefficient
Phi  = exp(F*Dt);               % state transition matrix (1x1)
Q    = Pinf*(1 - Phi^2);        % process noise covariance
SDQ  = sqrt(Q);                 % process noise standard deviation
R    = 1;                       % measurement noise covariance
SDR  = sqrt(R);                 % measurement noise standard deviation
%
% Create a simulated trajectory
%
x(1) = sqrt(Pinf)*randn;        % true initial value (random)
t1(1) = 0;                       % initial time
z(1) = H*x(1) + SDR*randn;      % measurement
for k=2:100,
    x(k) = Phi*x(k-1) + SDQ*randn;
    t1(k) = t1(k-1) + Dt;
    z(k) = H*x(k) + SDR*randn;
end;
%
% Forward filter pass
%
MSf0      = 0;                  % mean-squared innovations
xhatf0(1) = 0;                  % expected value for zero-mean process
Pf0(1)    = 1;                  % a priori variance
t(1)      = 0;                  % initial time
zhat      = H*xhatf0(1);        % expected measurement
Kbar      = Pf0(1)*H/(H*Pf0(1)*H+R);  % Kalman gain
xhatf0(2) = xhatf0(1) + Kbar*(z(1)-zhat);
MSf0      = MSf0 + (x(1)-xhatf0(2))^2;
Pf0(2)    = Pf0(1) - Kbar*H*Pf0(1);
t(2)      = 0;
j         = 2;
for k=2:100;
    j       = j + 1;
    t(j)    = t(j-1) + Dt;
    xhatf0(j) = Phi*xhatf0(j-1);    % predicted value
    Pf0(j)    = Phi*Pf0(j-1)*Phi + Q; % a priori variance
    zhat      = H*xhatf0(j);        % expected measurement
    Kbar      = Pf0(j)*H/(H*Pf0(j)*H+R);  % Kalman gain
    j         = j + 1;
    t(j)      = t(j-1);
    xhatf0(j) = xhatf0(j-1) + Kbar*(z(k)-zhat);
    MSf0      = MSf0 + (x(k)-xhatf0(j))^2;
    Pf0(j)    = Pf0(j-1) - Kbar*H*Pf0(j-1);
end;
RMSf0 = sqrt(MSf0/100);
figure;
subplot(3,1,1),
plot(t1,x,'k-',t,xhatf0,'k--',t,sqrt(Pf0),'k:',1,4.9,'w.',1,-4.9,'w.');
ylabel('(Forward)');
xlabel('Time');
title('Three-pass smoothing');
%
% Backward filter pass
%
MSf1        = 0;
xhatf1(200) = 0;                  % expected value for zero-mean process
Pf1(200)    = 1;                  % a priori variance
zhat        = H*xhatf1(200);        % expected measurement
Kbar        = Pf1(200)*H/(H*Pf1(200)*H+R);  % Kalman gain
xhatf1(199) = xhatf1(200) + Kbar*(z(100)-zhat);
MSf1        = MSf1 + (x(100)-xhatf1(199))^2;
Pf1(199)    = Pf1(200) - Kbar*H*Pf1(200);
j           = 199;
for k=99:-1:1;
    j       = j - 1;
    xhatf1(j) = xhatf1(j+1)/Phi;    % predicted value
    Pf1(j)    = Pf1(j+1)/Phi^2 + Q; % a priori variance
    zhat      = H*xhatf1(j);        % expected measurement
    Kbar      = Pf1(j)*H/(H*Pf1(j)*H+R);  % Kalman gain
    j         = j - 1;
    xhatf1(j) = xhatf1(j+1) + Kbar*(z(k)-zhat);
    MSf1      = MSf1 + (x(k)-xhatf1(j))^2;
    Pf1(j)    = Pf1(j+1) - Kbar*H*Pf1(j+1);
end;
RMSf1 = sqrt(MSf1/100);
subplot(3,1,2),
plot(t1,x,'k-',t,xhatf1,'k--',t,sqrt(Pf1),'k:',1,4.9,'w.',1,-4.9,'w.');
ylabel('(Backward)');
xlabel('Time');
%
% Smoothing pass
%
MSs1 = 0;
for k=1:100,
    Ps1(k) = 1/(1/Pf0(2*k)+1/Pf1(2*k));  % a post. forward, a pri. backward
    xs1(k) = (xhatf0(2*k)/Pf0(2*k)+xhatf1(2*k)/Pf1(2*k))*Ps1(k);
    MSs1   = MSs1 + (x(k) - xs1(k))^2;
    Ps2(k) = 1/(1/Pf0(2*k-1)+1/Pf1(2*k-1)); % pri. forward, post. backward
    xs2(k) = (xhatf0(2*k-1)/Pf0(2*k-1)+xhatf1(2*k-1)/Pf1(2*k-1))*Ps2(k);
end;
RMSs1 = sqrt(MSs1/100);
subplot(3,1,3),
plot(t1,x,'k-',t1,xs1,'k--',t1,sqrt(Ps1),'k:',1,4.9,'w.',1,-4.9,'w.');
ylabel('(Smoother)');
xlabel('Time');
legend('true','est.','\sigma');
disp(['RMS estimation error, forward filter pass  = ',num2str(RMSf0)]);
disp(['RMS estimation error, backward filter pass = ',num2str(RMSf1)]);
disp(['RMS estimation error, three-pass smoother  = ',num2str(RMSs1)]);
disp(' ');
disp('Results will differ with each Monte Carlo run.');
disp('(Try it.)');
