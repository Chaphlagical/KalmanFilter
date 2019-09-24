% 
% M. S. Grewal and A. P. Andrews,
% Kalman Filtering: Theory and Practice Using MATLAB, 4th Edition
% Wiley, 2014.
%
% Demonstration of probability conditioning on measurements of a random
% walk variate.
%
close all;
clear all;
disp('M. S. Grewal and A. P. Andrews,');
disp('Kalman Filtering: Theory and Practice Using MATLAB');
disp('4th Edition, Wiley ,2014.');
disp(' ');
disp('ProbCond.m:');
disp('Demonstration should help you visualize');
disp('how estimation influences uncertainty.');
disp(' ');
disp('It demonstrates the effect that the Kalman filter');
disp('has in terms of its conditioning of the probability');
disp('distribution of the state of a random walk variate.');
P = .4; x = sqrt(P)*randn(1); xhat = 0; H = 1; Q = .01; 
R = .04; A = zeros(101,101);nm = 0; mi = 21;
   for k=1:101,
   xt(k) = x;
   xh(k) = xhat;
   xu(k) = xhat+sqrt(P);
   xl(k) = xhat-sqrt(P);
      for m=1:101,
      dev  = (m-51)/25 - xhat;
      prob = exp(-dev^2/(2*P))/sqrt(2*pi*P);
      A(m,k) = prob;
      end;
   P = P + Q;
      if (mi*round(k/mi) == k),
      disp(['Measurement taken at t = ',num2str(k),' seconds.']);
      z = H*x + randn(1)*sqrt(R);
      K = P*H/(H^2*P+R);
      xhat = xhat + K*(z - H*xhat);
      P = P - K*H*P;
      nm = nm + 1;
      t0(nm) = k-1;
      z0(nm) = z;
      end;
   x = x + randn(1)*sqrt(Q);
   s(k) = (k-51)/25;
   t(k) = k-1;
   end;
clf;
surf(t,s,A);
zlabel('Probability Density');
xlabel('Time in Seconds');
ylabel('State Value');
title('Evolution of probability distribution of the state of a random walk');
disp('Plot shows evolution of inferred probability distribution');
disp('of the state of a random walk.');
disp(' ');
disp('Note the broadening of the probability distribution');
disp(['except when measurements are taken (every ',num2str(mi),' seconds).']);
disp(' ');
response = input('Press <ENTER> to continue.');
   if (response=='p'),                    % Allows printing the figure ...
   response = input('Print to what?','s');
     if (length(response)==7),
        if (response == 'printer'),print; % to the printer ...
        else print -deps response;end     % or to a file.
     else print -deps response;end;
   end;
figure;
plot(t,xt,t,xh,t,xu,t,xl);
   for n=1:nm,
   text(t0(n),z0(n),'z');
   end;
xlabel('Time in Seconds');
ylabel('State Value');
legend('true','est.','+1\sigma','-1\sigma');
title('True and estimated state of random walk.  "z" at measured values.');
disp('Plots show true and estimated (+/- 1 sigma) state values.');
disp('The true value should be within the +/- bounds');
disp(' more than half the time.  Is it?');
disp(' ');
disp('Measured values (with measurement noise)');
disp(' are indicated by "z"');
disp(' ');
disp('Why do the +/- 1-sigma bounds contract');
disp('whenever measurements are taken?');
disp('Why do they expand between measurements?');
disp(' ');
disp('NOTE: Results will be different each time you run this.');
disp('      Try it and see.');

