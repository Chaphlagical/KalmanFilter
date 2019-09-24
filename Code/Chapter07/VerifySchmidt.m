%
% Verify function Schmidt on pseudorandom matrices A(n,n), B(n,r)
% by computing C = Schmidt(A,B) and the relative error
% max(max(abs(A*A'+B*B'-C*C')))/max(max(abs(C*C')))
%
%   M. S. Grewal and A. P. Andrews,
%   Kalman Filtering: Theory and Practice Using MATLAB,
%   Wiley, 2008
%
clear all;
close all;
%
m1 = 0;
m2 = 0;
k  = 0;
q  = 12;
disp('Verify function Schmidt on pseudorandom matrices A(n,n), B(n,r)');
disp(['  for 1 <= n <= ',num2str(q),' and 1 <= r <= ',num2str(q)]);
disp('by computing C = Schmidt(A,B), the relative solution error');
disp('   max(max(abs(A*A''+B*B''-C*C'')))/max(max(abs(C*C'')))');
disp('and the relative triangularization error.');
for n=1:q,
    for r=1:q,
        A = randn(n);
        B = randn(n,r);
        C = Schmidt(A,B);
        D = A*A'+B*B';
        E = C*C';
        d = max(max(abs(D-E)))/max(max(abs(E)));
        upperss  = 0;
        upperno  = 0;
        lowerss  = 0;
        lowerno  = 0;
        if n>1
            for j=1:n,
                for i=1:n,
                    if i>j
                        lowerss = lowerss + C(i,j)^2;
                        lowerno = lowerno + 1;
                    else
                        upperno  = upperno + 1;
                        upperss  = upperss + C(i,j)^2;
                    end;
                end;
            end;
        perfmeas = sqrt(lowerss/lowerno)/sqrt(upperss/upperno);
        disp(['n = ',num2str(n),', r = ',num2str(r),', rel.sol.err. = ',num2str(d),' rel.tri.err = ',num2str(perfmeas)]);
        m2 = max(m2,perfmeas);
        end;
        m1 = max(m1,d);
        k  = k+1;
    end;
end;
disp(['Maximum relative solution error on ',num2str(k),' pseudorandom samples = ',num2str(m1)]);
disp(['Maximum triangularization error on ',num2str(k),' pseudorandom samples = ',num2str(m2)]);
