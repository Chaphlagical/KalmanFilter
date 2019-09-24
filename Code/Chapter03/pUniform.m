function px = pUniform(x,mean,variance)
%
% returns the probability density at x of a uniform distribution
% with given meanand variance
%
mu    = mean;
sigma = sqrt(variance);
if x < mu - sqrt(3)*sigma,
    px  = 0;
elseif x > mu + sqrt(3)*sigma,
    px  = 0;
else
    px  = 1/2/sqrt(3)/sigma;
end;