function H = dhdxarctan(x,d);
%
% Arctangent measurement sensitivity function
%
H = d/(d^2 + x^2);