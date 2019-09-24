%
% Example 5.8
%
close all;
clear all;
clf;
disp('Covariance analysis of radar tracking problem,');
disp('given as Example 5.8 in');
disp('M. S. Grewal and A. P. Andrews,');
disp('Kalman Filtering: Theory and Practice Using MATLAB');
disp('4th Edition, Wiley, 2014.');
disp(' ');
disp('Plots histories of six mean squared state');
disp('uncertainties and six magnitudes of Kalman gains');
disp('for intersample intervals of 5, 10 and 15 seconds.');
disp(' ');
disp('Six state variables:');
disp('  1. Range to object being tracked.');
disp('  2. Range rate of object being tracked.');
disp('  3. Object range maneuvering noise (pseudo state).');
disp('  4. Bearing to object being tracked.');
disp('  5. Bearing rate of object being tracked.');
disp('  6. Object bearing maneuvering noise (pseudo state).');
disp('Pseudo states are used for modeling correlated noise.');
%
sigma1sq = (103/3)^2;
sigma2sq = 1.3E-8;
sigmarsq = (1000)^2;
sigmatsq = (.017)^2;
rho      = 0.5;
%
% State Transition matrix (the part not depending on T)
%
Phi      = eye(6);
Phi(2,3) = 1;
Phi(5,6) = 1;
Phi(3,3) = rho;
Phi(6,6) = rho;
Q        = zeros(6);
Q(3,3)   = sigma1sq;
Q(6,6)   = sigma2sq;
R        = zeros(2);
R(1,1)   = sigmarsq;
R(2,2)   = sigmatsq;
H        = zeros(2,6);
H(1,1)   = 1;
H(2,4)   = 1;
%
% arrays for saving data to be plotted
%
t        = zeros(3,32); % time (for 3 plots)
rcov     = zeros(3,32); % Range covariance
rrcov    = zeros(3,32); % Range Rate covariance
bcov     = zeros(3,32); % Bearing covariance
brcov    = zeros(3,32); % Bearing Rate covariance
rrncov   = zeros(3,32); % Range Rate Noise covariance
brncov   = zeros(3,32); % Bearing Rate Noise covariance
rkg      = zeros(3,32); % Range Kalman gain
rrkg     = zeros(3,32); % Range Rate Kalman gain
bkg      = zeros(3,32); % Bearing Kalman gain
brkg     = zeros(3,32); % Bearing Rate Kalman gain
rrnkg    = zeros(3,32); % Range Rate Noise Kalman gain
brnkg    = zeros(3,32); % Bearing Rate Noise Kalman gain
N=0;
   for T = 5:5:15,
   N=N+1;
   disp(['Simulating tracking at ',num2str(T),' second intervals.']);
   Phi(1,2) = T;
   Phi(4,5) = T;
   P        = zeros(6);
   P(1,1)   = sigmarsq;
   P(1,2)   = sigmarsq/T;
   P(2,1)   = P(1,2);
   P(2,2)   = 2*sigmarsq/T^2 + sigma1sq;
   P(3,3)   = sigma1sq;
   P(4,4)   = sigmatsq;
   P(4,5)   = sigmatsq/T;
   P(5,4)   = P(4,5);
   P(5,5)   = 2*sigmatsq/T^2 + sigma2sq;
   P(6,6)   = sigma2sq;
      for cycle=0:15,
%
% Save a priori values
%
      prior      = 2*cycle+1;
      t(N,prior) = T*cycle;
      K          = P*H'/(H*P*H'+R);
      rcov(N,prior)   = P(1,1); % Range covariance
      rrcov(N,prior)  = P(2,2); % Range Rate covariance
      bcov(N,prior)   = P(4,4); % Bearing covariance
      brcov(N,prior)  = P(5,5); % Bearing Rate covariance
      rrncov(N,prior) = P(3,3); % Range Rate Noise covariance
      brncov(N,prior) = P(6,6); % Bearing Rate Noise covariance
      rkg(N,prior)    = sqrt(K(1,1)^2+K(1,2)^2); % Range Kalman gain
      rrkg(N,prior)   = sqrt(K(2,1)^2+K(2,2)^2); % Range Rate Kalman gain
      bkg(N,prior)    = sqrt(K(4,1)^2+K(4,2)^2); % Bearing Kalman gain
      brkg(N,prior)   = sqrt(K(5,1)^2+K(5,2)^2); % Bearing Rate Kalman gain
      rrnkg(N,prior)  = sqrt(K(3,1)^2+K(3,2)^2); % Range Rate Noise Kalman gain
      brnkg(N,prior)  = sqrt(K(6,1)^2+K(6,2)^2); % Bearing Rate Noise Kalman gain
%
% Save a posteriori values
%
      P         = P - K*H*P;
      P         = .5*(P+P');
      post      = prior + 1;
      t(N,post) = T*cycle;
      rcov(N,post)   = P(1,1); % Range covariance
      rrcov(N,post)  = P(2,2); % Range Rate covariance
      bcov(N,post)   = P(4,4); % Bearing covariance
      brcov(N,post)  = P(5,5); % Bearing Rate covariance
      rrncov(N,post) = P(3,3); % Range Rate Noise covariance
      brncov(N,post) = P(6,6); % Bearing Rate Noise covariance
      rkg(N,post)    = sqrt(K(1,1)^2+K(1,2)^2); % Range Kalman gain
      rrkg(N,post)   = sqrt(K(2,1)^2+K(2,2)^2); % Range Rate Kalman gain
      bkg(N,post)    = sqrt(K(4,1)^2+K(4,2)^2); % Bearing Kalman gain
      brkg(N,post)   = sqrt(K(5,1)^2+K(5,2)^2); % Bearing Rate Kalman gain
      rrnkg(N,post)  = sqrt(K(3,1)^2+K(3,2)^2); % Range Rate Noise Kalman gain
      brnkg(N,post)  = sqrt(K(6,1)^2+K(6,2)^2); % Bearing Rate Noise Kalman gain
      P              = Phi*P*Phi' + Q;
      end;
   end;
%figure;
subplot(2,2,1),plot(t(1,:),rcov(1,:),t(2,:),rcov(2,:),t(3,:),rcov(3,:));xlabel('Time (sec)');ylabel('Range Cov');
subplot(2,2,2),plot(t(1,:),rrcov(1,:),t(2,:),rrcov(2,:),t(3,:),rrcov(3,:));xlabel('Time (sec)');ylabel('Range Rate Cov');
subplot(2,2,3),plot(t(1,:),bcov(1,:),t(2,:),bcov(2,:),t(3,:),bcov(3,:));xlabel('Time (sec)');ylabel('Bear. Cov');
subplot(2,2,4),plot(t(1,:),brcov(1,:),t(2,:),brcov(2,:),t(3,:),brcov(3,:));xlabel('Time (sec)');ylabel('Bear. Rate Cov');
disp('Figure 5.10');
disp('Press <ENTER> for more plots.');
pause
figure;
subplot(2,2,1),plot(t(1,:),rrncov(1,:),t(2,:),rrncov(2,:),t(3,:),rrncov(3,:));xlabel('Time (sec)');ylabel('RRateNoise Cov');
subplot(2,2,2),plot(t(1,:),brncov(1,:),t(2,:),brncov(2,:),t(3,:),brncov(3,:));xlabel('Time (sec)');ylabel('BRateNoise Cov');
subplot(2,2,3),plot(t(1,:),rkg(1,:),t(2,:),rkg(2,:),t(3,:),rkg(3,:));xlabel('Time (sec)');ylabel('Range Gain');
subplot(2,2,4),plot(t(1,:),rrkg(1,:),t(2,:),rrkg(2,:),t(3,:),rrkg(3,:));xlabel('Time (sec)');ylabel('Range Rate GAin');
disp('Figure 5.11');
disp('Press <ENTER> for more plots.');
pause
figure;
subplot(2,2,1),plot(t(1,:),bkg(1,:),t(2,:),bkg(2,:),t(3,:),bkg(3,:));xlabel('Time (sec)');ylabel('Bear. Gain');
subplot(2,2,2),plot(t(1,:),brkg(1,:),t(2,:),brkg(2,:),t(3,:),brkg(3,:));xlabel('Time (sec)');ylabel('Bear Rate Gain');
subplot(2,2,3),plot(t(1,:),rrnkg(1,:),t(2,:),rrnkg(2,:),t(3,:),rrnkg(3,:));xlabel('Time (sec)');ylabel('RRateNoiseGain');
subplot(2,2,4),plot(t(1,:),brnkg(1,:),t(2,:),brnkg(2,:),t(3,:),brnkg(3,:));xlabel('Time (sec)');ylabel('BRateNoiseGain');
disp('Figure 5.12');
disp('DONE');