%
disp('GREWAL AND ANDREWS,');
disp('KALMAN FILTERING: THEORY AND PRACTICE USING MATLAB, 4TH EDITION,');
disp('WILEY, 2014.');
disp('CHAPTER 3');
disp('--------------------------------------------------------------------');
disp('meansqesterr.m');
disp('--------------------------------------------------------------------');
disp('Computes the mean-squared estimation error as a function of the estimate');
disp('for five probability distributions: Gaussian, chi-squared (2dof), Laplace,');
disp('and log-normal.  The mean is 2 (two) and the variance is 4 (four) in all');
disp('cases, so all probability density functions will have the same mean and ');
disp('variance.');
disp('--------------------------------------------------------------------');
disp('Be prepared: RUN TIME IS LIKELY TO BE A FEW MINUTES.');
disp('--------------------------------------------------------------------');
%
% NOTE1:The lognormal distribution has what is called a "heavy tail,"
%       meaning that it decays sub-exponentially as x -> +infinity.
%
%       As a consequence, it is necessary to include a domain that is 
%       orders of magnitude greater than its standard deviation in order
%       to get reasonably accurate results.
%
%       In the example below with variance equal to 4 (standard deviation
%       equal to 2) the range of the variate used in integration is from
%       16*std below the mean to 31*std above the mean, because the fat
%       tail is only on the positive side.  There is no tail on the
%       negative side.
%
% NOTE2: The chi-squared distribution has only one parameter, called its
%        "degrees of freedom,"  which is also its mean.  Therefore, using
%        a mean value of all distributions equal to 2 (two) results in a
%        mean of 2 (two) and a variance of 4 (four) for the chi-squared 
%        distribution.  The mean was then set to 2 and the variance to 4 
%        for all distributions.
%
% Functions performed:
%
%   1. Use the associated functions
%       1.1 pNormal for normal (or Gaussian) probability densities.
%       1.2 pLaplace for LaPlace probability densities.
%       1.3 pLogNormal for lognormal probability densities.
%       1.4 pChiSquared chi-squared probability densities.
%       1.5 pUniform for uniform probability densities.
%      to calculate the probability density function values over the range
%      from -16 standard deviations to +31 standard deviations from the
%      mean.
%   2. Plot the functions over the range: mean +/- 3 standard deviations
%   3. As a test of the functions and parameters used, calculate and
%      display the computed zeroth and first raw moments of the computed 
%      distributions, and their 2nd central moments (variances).
%   4. Calculate and plot the computed values of the mean-squared
%      estimation error as a function of the value of the estimate.  These
%      plots should show a minimum at the mean of each distribution.
%
close all;
clear all;
tic;
mean     = 2;
variance = 2*mean;
%
% 1. Calculate the probability distributions at intervals dx
%    (Because the number of values of x used is 94001, this can take a
%    while.)
%
disp('(Probabilities being computed.)');
disp('--------------------------------------------------------------------');
k        = 0;
dx       = 0.001;
for x=-30:dx:64,
    k       = k + 1;
    xk(k)   = x;
    pNo(k)  = pNormal(x,mean,variance);
    pLp(k)  = pLaplace(x,mean,variance);
    pLN(k)  = pLogNormal(x,mean,variance);
    pXs(k)  = pChiSquared(x,mean);
    pUn(k)  = pUniform(x,mean,variance);
end;
%
disp('(Probabilities computed, being plotted.)');
disp('--------------------------------------------------------------------');
toc;
%
%
% 2. Plot distribution functions from x=-6 to x=10: (-4*std to + 4*std)
%
k1  = 24/dx+1;
k2  = 40/dx+1;
%
figure;
plot(xk(k1:k2),pUn(k1:k2),'k-.',xk(k1:k2),pNo(k1:k2),'k-',xk(k1:k2),pLp(k1:k2),'k--',xk(k1:k2),pLN(k1:k2),'k-',xk(k1:k2),pXs(k1:k2),'k:',[mean-sqrt(variance),mean,mean+sqrt(variance)],[.6,.6,.6],'k-+','LineWidth',2);
xlabel('x','FontSize',18);
ylabel('p(x)','FontSize',18);
text(mean,.63,'MEAN \pm 1{\times}STD','HorizontalAlignment','center','FontSize',12);
title('FIVE PROBABILITY DENSITIES WITH THE SAME MEAN AND VARIANCE','FontSize',12);
legend('Uniform','Normal','Laplace','logNormal','ChiSq2dof');
disp('Probabilities plotted');
toc;
%
% 3. Calculate and display the moments of order 
%       0, should equal 1 in all cases
%       1, should equal the mean, 2 in all cases
%       2nd central moment, should equal the variance, 4 in all cases
%
%    (Because the number of values of x used is 94001, this can take a
%    while.)
%
kmax    = length(x);
xk      = xk';
pNo     = pNo';
pLp     = pLp';
pLN     = pLN';
pXs     = pXs';
pUn     = pUn';
xD      = xk - 2*ones(kmax,1);  % Central difference
OneNo   = dx*sum(pNo);
OneLp   = dx*sum(pLp);
OneLN   = dx*sum(pLN);
OneXs   = dx*sum(pXs);
OneUn   = dx*sum(pUn);
meanNo  = dx*xk'*pNo;
meanLp  = dx*xk'*pLp;
meanLN  = dx*xk'*pLN;
meanXs  = dx*xk'*pXs;
meanUn  = dx*xk'*pUn;
varNo   = dx*pNo'*xD.^2;
varLp   = dx*pLp'*xD.^2;
varLN   = dx*pLN'*xD.^2;
varXs   = dx*pXs'*xD.^2;
varUn   = dx*pUn'*xD.^2;
disp('Moments computed');
disp('--------------------------------------------------------------------');
toc;
disp(' ');
disp('Check distribution statistical parameters by rectangular numerical integration:');
disp('--------------------------------------------------------------------');
disp('DISTRIBUTION    _____________ MOMENTS ___________');
disp('                    0th       1st     2nd Central');
disp(['Uniform             ',num2str(OneUn,3),'         ',num2str(meanUn,3),'       ',num2str(varUn,3)]);
disp(['Normal              ',num2str(OneNo,3),'         ',num2str(meanNo,3),'       ',num2str(varNo,3)]);
disp(['Laplace             ',num2str(OneLp,3),'         ',num2str(meanLp,3),'       ',num2str(varLp,3)]);
disp(['ChiSquared w/ 2DOF) ',num2str(OneXs,3),'         ',num2str(meanXs,3),'       ',num2str(varXs,3)]);
disp(['LogNormal           ',num2str(OneLN,3),'         ',num2str(meanLN,3),'       ',num2str(varLN,3)]);
disp('--------------------------------------------------------------------');
%
%   4. Calculate and plot the mean-squared estimation error as a function
%      of the value of the estimate.  
%
xhat0  = ones(kmax,1);
%
k      = 0;
disp('Calculating squared errors.');
disp('--------------------------------------------------------------------');
for i=k1:k2,
    k   = k + 1;
    if mod(k,1000)==1, 
        disp(['k = ',num2str(k),' of ',num2str(k2-k1+1)]);
        toc;
    end;
    xhat        = xk(i)*xhat0;
    SqError     = (xhat - xk).^2;
    SqErNo(k)   = dx*pNo'*SqError;
    SqErLp(k)   = dx*pLp'*SqError;
    SqErLN(k)   = dx*pLN'*SqError;
    SqErXs(k)   = dx*pXs'*SqError;
    SqErUn(k)   = dx*pUn'*SqError;
end;
x   = xk(k1:k2);
disp('Squared estimation errors computed');
disp('--------------------------------------------------------------------');
toc;
figure;
semilogy(x,SqErUn,'k-.',x,SqErNo,'k-',x,SqErLp,'k--',x,SqErLN,'k-',x,SqErXs,'k:','LineWidth',2);
xlabel('Estimate Value','FontSize',18);
ylabel('Mean-Squared Estimation Error','FontSize',18);
legend('Uniform','Normal','Laplace','logNormal','ChiSq2dof');
text(2,2,'(Mean-square estimation error should be minimum (4) at mean (2).)','HorizontalAlignment','center');
disp('Squared estimation errors plotted');
disp('--------------------------------------------------------------------');
toc;
figure;
MinSqEr     = min(min([SqErLN;SqErLp;SqErNo;SqErUn;SqErXs]));
plot(x,(SqErUn-SqErLN)/MinSqEr,'k-',x,(SqErLp-SqErLN)/MinSqEr,'k--',x,(SqErXs-SqErLN)/MinSqEr,'k-.',x,(SqErNo-SqErLN)/MinSqEr,'k:','LineWidth',2);
ylabel('Relative Squared Error','FontSize',18);
xlabel('Estimate Value','FontSize',18);
legend('SqErUn-SqErLN','SqErLp-SqErLN','SqErXs-SqErLN','SqErNo-SqErLN');
disp('Squared estimation error differences plotted');
disp('--------------------------------------------------------------------');
toc;
disp('Note that their values at \hat{x} = 2 are not exactly zero,');
disp('which is due to the limited accuracy of numerical integration,');
disp('---especially for the "high-tailed" lognormal distribution.');
disp('--------------------------------------------------------------------');
%