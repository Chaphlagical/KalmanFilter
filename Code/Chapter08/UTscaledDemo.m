%
% Scaled unscented transform with eigenvector Cholesky factors
%
% 
% M. S. Grewal & A. P. Andrews
% Kalman Filtering Theory and Practice Using MATLAB
% 4th Edition, Wiley, 2014.
% 
close all;
clear all;
%
% Initial covariance
%
T               = 1/2*[1,1,1,1;1,1,-1,-1;1,-1,1,-1;1,-1,-1,1]; % orthogonal
P               = T*[4,0,0,0;0,36,0,0;0,0,100,0;0,0,0,196]*T';
%
% Eigenvector Choleksy factor
%
[U,Lambda,V]    = svd(P);
CEV             = U*sqrt(Lambda);
%
% Initial mean
%
xhat = [0;0;0;0];
%
% Theoretical mean from NonLin transformation of (xhat(1:3),P(1:3,1:3))
%
Tzhat    = zeros(3,1);
Tzhat(1) = P(2,3)*(1 - P(2,3)^2 / P(2,2) / P(3,3))^(-3/2);
Tzhat(2) = P(3,1)*(1 - P(3,1)^2 / P(3,3) / P(1,1))^(-3/2);
Tzhat(3) = P(1,2)*(1 - P(1,2)^2 / P(1,1) / P(2,2))^(-3/2);
%
% Extendeded Kalman filter transformation of mean and variance
%
Lzhat = NonLin(xhat(1:3));
H     = [0,xhat(3),xhat(2),0;xhat(3),0,xhat(1),0;xhat(2),xhat(1),0,0];
LPzz  = H*P*H';
%
% Unscented Transform solutions with kappa = 0...2
%
kappa = 0;
alpha = 1;
beta  = 2;
[sUTzhatEV,sUTPzzEV]    = UTscaled(xhat(1:3),CEV(1:3,1:3),@NonLin,kappa,alpha,beta);
k = 0;
yhat = [];
for kappa=0:.1:2,
    k = k + 1;
    [sUTzhatEV,sUTPzzEV]    = UTscaled(xhat(1:3),CEV(1:3,1:3),@NonLin,kappa,alpha,beta);
    Kappa(k) = kappa;
    yhat = [yhat,sUTzhatEV];
end;
figure;
plot(Kappa,yhat(1,:),'k-',Kappa,yhat(2,:),'k--',Kappa,yhat(3,:),'k:');
hold on;
plot([0,2],[Tzhat(1),Tzhat(1)],'b-',[0,2],[Tzhat(2),Tzhat(2)],'b--',[0,2],[Tzhat(3),Tzhat(3)],'b:');
hold off;
legend('z_1','z_2','z_3');
xlabel('\kappa');
title('UT solution vs \kappa (true values in blue)');
ylabel('z-hat');
%
% Vary alpha from 0.05 to 1
%
kappa = 1;
beta  = 2;
k = 0;
yhat = [];
for alpha=0.05:0.05:1,
    k = k + 1;
    [sUTzhatEV,sUTPzzEV]    = UTscaled(xhat(1:3),CEV(1:3,1:3),@NonLin,kappa,alpha,beta);
    Alpha(k) = alpha;
    yhat = [yhat,sUTzhatEV];
end;
figure;
plot(Alpha,yhat(1,:),'k-',Alpha,yhat(2,:),'k--',Alpha,yhat(3,:),'k:');
hold on;
plot([0,1],[Tzhat(1),Tzhat(1)],'b-',[0,1],[Tzhat(2),Tzhat(2)],'b--',[0,1],[Tzhat(3),Tzhat(3)],'b:');
hold off;
legend('z_1','z_2','z_3');
xlabel('\alpha');
title('UT solution vs \alpha (true values in blue)');
ylabel('z-hat');
 
