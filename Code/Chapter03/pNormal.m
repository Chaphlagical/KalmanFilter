function px = pNormal(x,mean,variance)
%
% returns the probability density at x of the normal distribution
% with given meanand variance
%
mu    = mean;
sigma = sqrt(variance);
if length(x) == 1,
    px    = exp(-(x-mu)^2/(2*sigma^2))/sqrt(2*pi*sigma^2);
else
    for i=1:length(x),
        px(i)    = exp(-(x(i)-mu)^2/(2*sigma^2))/sqrt(2*pi*sigma^2);
    end;
end;