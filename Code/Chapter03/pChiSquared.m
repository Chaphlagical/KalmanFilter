function px = pChiSquared(x,mean)
%
% returns the probability density at x of the chi-squared distribution
% with given mean.
%
% For chi-squared distributions, variance is always twice the mean.
%
k   = mean;  % degrees of freedom
if x>=0,
    px  = x^(k/2-1)*exp(-x/2)/2^(k/2)/gamma(k/2);
else
    px  = 0;
end;
