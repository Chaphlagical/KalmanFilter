function px = pLogNormal(x,mean,variance)
%
% returns the probability density at x of the lognormal distribution
% with given mean and variance
%
% NOTE: The lognormal distribution has what is called a "heavy tail,"
%       meaning that it decays sub-exponentially as x -> infinity.
%
%       As a consequence, it is necessary to include a domain that is 
%       orders of magnitude greater than its standard deviation in order
%       to get reasonably accurate results.
%
sigma = sqrt(log(1+variance/mean^2));
mu    = log(mean^2/sqrt(variance+mean^2));
if x>0,
    px    = exp(-(log(x)-mu)^2/2/sigma^2)/x/sigma/sqrt(2*pi);
else
    px  = 0;
end;