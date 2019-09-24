function p = pYArctanaX(y,a)
%
% Computes the probability density function of
%
%      Y = arctan(a X)
%
% where X has a zero-mean unit-normal univariate distribution.
%
if abs(y) >= pi/2
    p = 0;
else
    p = (1 + tan(y)^2)/a*sqrt(2)*pi^(-1/2)*exp(-(tan(y)/a)^2/2)/2;
end;