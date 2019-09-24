disp('M. S. Grewal and A. P. Andrews,');
disp('Kalman Filtering: Theory and Practice using MATLAB, 4th Edition,');
disp('Wiley, 2014.');
disp(' ');
disp('Calculates and plots RMS estimation uncertainty of fixed-point smoother');
disp('for a scalar linear time-invariant random process in continuous time,');
disp('with specified model parameters.');
disp(' ');
%
close all;
clear all;
F  = 1/100;  % dynamic coefficient of system dynamic model
Q  = 1/100;  % variance of dynamic process noise
H  = 1;      % measurement sensitivity
R  = 0.01;   % variance of measurement noise
t0 = 0;      % starting time
tF = 100;    % fixed time at which system state is to be estimated
te = 300;    % ending time
%
k = 0;
for t=t0+1:te,
    k     = k + 1;
    ts(k) = t;
    if t<tF
        tf(k)  = t;
        rho    = sqrt(R*(F^2*R + Q*H^2));
        Eminus = exp(-rho/R*(t-t0));
        Eplus  = exp(rho/R*(t-t0));
        t2     = R*F;
        t7     = H^2;
        Pf(k)  = (Eplus*rho-t2*Eminus+t2*Eplus+Eminus*rho)/t7/(Eplus-Eminus);
        t4     = exp(2*F*(tF - t));
        t5     = Pf(k);
        Ps(k)  = t4*t5+Q*(t4-1)/F/2;
        Ps0    = Ps(k);
        %
    elseif t==tF
        rho    = sqrt(R*(F^2*R + Q*H^2));
        Eminus = exp(-rho/R*(t-t0));
        Eplus  = exp(rho/R*(t-t0));
        t2     = R*F;
        t7     = H^2;
        Ps(k)  = (Eplus*rho-t2*Eminus+t2*Eplus+Eminus*rho)/t7/(Eplus-Eminus);
        tf(k)  = t;
        Pf(k)  = Ps(k);
        Ps0    = Ps(k);
    else
        tb(k)  = t;
        rho    = sqrt(R*(F^2*R + Q*H^2));
        Eminus = exp(-rho/R*(t-tF));
        Eplus  = exp(rho/R*(t-tF));
        t2     = R*F;
        t7     = H^2;
        Pb(k)  = (Eminus+Eplus)*rho/t7/(Eplus-Eminus)-t2/t7;
        Ps(k)  = Ps0*Pb(k)/(Ps0+Pb(k));
    end;
end;
semilogy(ts,sqrt(Ps),'k-',tf,sqrt(Pf),'k--',tb,sqrt(Pb),'k:','LineWidth',1.5);
legend('FP smooth','Fwd filt','Bkwd filt');
xlabel('Time [sec]','FontSize',15);
ylabel('RMS smoother uncertainty','FontSize',15);
title('Scalar linear time-invariant continuous fixed-point smoother model');
disp(['F  = ',num2str(F),', (dynamic coefficient)']);
disp(['Q  = ',num2str(Q),', (variance of dynamic process noise)']);
disp(['H  = ',num2str(H),', (measurement sensitivity)']);
disp(['R  = ',num2str(R),', (variance of measurement noise)']);
disp(['t0 = ',num2str(t0),', (starting time)']);
disp(['tF = ',num2str(tF),', (fixed time at which system state is to be estimated)']);
disp(['te = ',num2str(te),', (ending time)']);
