close all;
clear all;
%
disp('----------------------------------------------------------------');
disp('GREWAL AND ANDREWS,');
disp('KALMAN FILTERING: THEORY AND PRACTICE USING MATLAB, 4TH EDITION,');
disp('WILEY, 2014.');
disp('----------------------------------------------------------------');
disp('CHAPTER 2');
%
% Demonstrating matrix exponential function
%
disp('----------------------------------------------------------------');
disp('Exponentials of randomly generated matrices N = expm(M)');
disp('----------------------------------------------------------------');
for n =1:6,
   M = randn(n)
   N = expm(M)
end;
disp('----------------------------------------------------------------');
disp('Exponential of scaled shift matrix');
disp('----------------------------------------------------------------');
n = 6;
s = 1/10;
M = zeros(n);
for k=2:n,
   M(k-1,k) = s;
end;
M
expm(M)
disp('----------------------------------------------------------------');
disp('Except for shift matrix, results should be different when expm1 is re-run.');
disp('----------------------------------------------------------------');
