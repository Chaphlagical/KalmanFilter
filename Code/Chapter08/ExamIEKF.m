%
% 
disp('M. S. Grewal & A. P. Andrews');
disp('Kalman Filtering Theory and Practice Using MATLAB');
disp('4th Edition, Wiley, 2014.');
disp('Chapter 8');
disp('Iterated Extended Kalman Filter (IEKF)');
disp('Application to nonlinear (arctangent) sensor');
disp(' ');
%
% Should demonstrate a reduction in errors due to nonlinearity
% by one or two orders of magnitude, using pseudorandom examples.
%
close all;
clear all;
%
d      = 1;
R      = 1e-6;   % 1 milliradian RMS sensor noise
RMS    = 1/2;    % RMS initial uncertainty
P      = RMS^2;
k      = 0;
figure;
xtrue  = 1;
semilogy([0,0],[.00000001,10000],'k-',[0,5],[.00000001,.00000001],'k-');
hold on;
for sample=1:100,
    k       = k + 1;
    z       = harctan(xtrue,1) + sqrt(R)*randn;
    xhat    = xtrue + RMS*randn; % 1-sigma initial error
    xnew    = xhat;
    Y       = [xhat-xtrue];
    X       = [0];
    for m=1:5; % Iteration counter
        H       = dhdxarctan(xnew,d);
        PHT     = P*H';
        K       = PHT/(H*PHT+R);
        dz      = z - harctan(xnew,d);
        xold    = xnew;
        xnew    = xhat + K*(dz - H*(xhat - xold));
        dx      = xnew - xold;
        dif     = dx'/P*dx;
        Y       = [Y,xnew-xtrue];
        X       = [X,m];
    end;
    semilogy(X,abs(Y),'k:');
end;
hold off;
xlabel('ITERATION NUMBER');
ylabel('ABSOLUTE ESTIMATION ERROR');
text(2.5,3e3,'PLOTS = 100 MONTE CARLO SIMULATIONS','HorizontalAlignment','center');
text(2.5,6e2,'USING ITERATED EXTENDED KALMAN FILTER','HorizontalAlignment','center');
text(2.5,1.2e2,'ON ARCTANGENT MEASUREMENT PROBLEM','HorizontalAlignment','center');
disp('---------------------------------------------------------------------');
disp('Re-running should give differrent Monte Carlo results.');
disp('Try it.');