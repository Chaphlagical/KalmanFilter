%
% PERFORMANCE EVALUATION OF FIXED INTERVAL SMOOTHER ON FINITE INTERVAL
%
disp('M. S. Grewal and A. P. Andrews,');
disp('Kalman Filtering: Theory and Practice using MATLAB, 4th Edition,');
disp('Wiley, 2014.');
disp(' ');
disp('Solves Riccati equations for a fixed-interval smoother');
disp('for the forward pass, backward pass, and combined covariances.');
disp(' ');
disp('See how it handles end effects (script), and the results (plotted).');
disp(' ');
%
close all;
clear all;
%
R = .01;
Q = .01;
H = 1;
F = .3;
%
tstart = 0;
tend   = 100;
%
k = 0;
for t=tstart:.5:tend,
    k  = k + 1;
    tk(k) = t;
    t1 = F ^ 2;
    t3 = H ^ 2;
    t7 = sqrt(R * (t1 * R + t3 * Q));
    tau = R / t7;
    if t==tstart
        Pf(k) = Inf;
    else
        t10 = (t-tstart)/tau;
        t11 = exp(t10);
        t12 = exp(-t10);
        t16 = -t12 + t11;
        Pf(k) = t7 * (t11 + t12) / t3 / t16 + F * R / t3;
    end
    %
    if t==tend
        Pb(k) = Inf;
    else
        t10 = (tend-t)/tau;
        t11 = exp(t10);
        t12 = exp(-t10);
        t15 = 0.1e1 / t3;
        Pb(k) = t7 * (t11 + t12) * t15 / (t11 - t12) - F * R * t15;
    end
    if t==tstart
        Ps(k) = Pb(k);
    elseif t==tend
        Ps(k) = Pf(k);
    else
        Ps(k) = Pf(k)*Pb(k)/(Pf(k)+Pb(k));
    end
end
plot(tk,Ps,'k-',tk,Pf,'k--',tk,Pb,'k:');
xlabel('FIXED TIME INTERVAL');
ylabel('VARIANCE OF ESTIMATION ERROR');
title(['F = ',num2str(F),', Q = ',num2str(Q),', H = ',num2str(H),', R = ',num2str(R)]);
midtime = (tstart+tend)/2;
midk    = round(length(tk)/2);
text(midtime,Ps(midk),'FIXED INTERVAL SMOOTHER','VerticalAlignment','bottom','HorizontalAlignment','center');
text(midtime,Pf(midk),'FORWARD FILTER','VerticalAlignment','bottom','HorizontalAlignment','center');
text(midtime,Pb(midk),'BACKWARD FILTER','VerticalAlignment','bottom','HorizontalAlignment','center');
disp('You can alter the values used for F, Q, H, and R, if you like,');
disp('and see the  results.');


    
