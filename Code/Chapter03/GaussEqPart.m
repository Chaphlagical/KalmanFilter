%   M. S. Grewal and A. P. Andrews
%   Kalman Filtering Theory and Practice Using MATLAB, 4th Edition,
%   Wiley, 2014
%   Chapter 3
%
function medians = GaussEqPart(n)
%
% Unit normal "democratic" partitioning.
%
% Parititions the domain of the zero-mean univariate Gaussian probability
% function into n segments, each of which has probability measure 1/n, and
% returns a vector of the medians of those segments.
%
% INPUT: n, the number of partitions
% 
% OUTPUT: x, vector of the medians of the respective segments
%
k   = 0;
for P = 1/2/n:1/n:1-1/2/n,
    k           = k+1;
    medians(k) = sqrt(2)*erfinv(2*P-1);
end;