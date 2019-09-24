function output = NonLin(input);
%
% Nonlinear vector function with output(j) =
%                                    input(mod(j-2,L)+1)*input(mod(j,L)+1);
% For example,
%   output(1) = input(rows)*input(2) 
%   output(2) = input(1)*input(3)
%
[L,M] = size(input);
for k=1:M,
    for j=1:L,
        output(j,k) = input(mod(j-2,L)+1,k)*input(mod(j,L)+1,k);
    end;
end;
return;