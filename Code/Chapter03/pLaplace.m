function px = pLaplace(x,mean,variance)
%
% returns the probability density at x of the Laplace distribution
% with given mean and variance
%
mu    = mean;
b     = sqrt(variance/2);
px    = exp(-abs(x-mu)/b)/(2*b);