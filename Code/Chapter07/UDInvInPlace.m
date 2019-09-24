%
disp('Grewal & Andrews,');
disp('Kalman Filtering Theory and Practice Using MATLAB, 4th Edition,');
disp('Wiley, 2014.');
disp(' ');
disp('Test of algorithms for UD decomposition and inversion of symmetric');
disp('positive definite matrix in place, overwriting the input with the output.');
disp(' ');
%
disp('Test matrices of dimension 1 to 12 are generated using random number');
disp('generator randn, with the number of samples equal to the dimension');
disp('of the matrix to be generated.');
disp('(This can result in deficient pseudorank.)');
disp(' ');
%
for m=1:12,
    M = zeros(m);
    X = zeros(m);
    for k=1:m,
        v = randn(m,1);
        M = M + v*v';
    end;
    UD = M;
    UD = UD_decomp(UD);
    for i=1:m,
        for j=i:m,
            if i==j
                s = UD(i,i);
            else
                s = UD(i,j)*UD(j,j);
            end;
            for k=j+1:m,
                s = s + UD(i,k)*UD(k,k)*UD(j,k);
            end;
            X(i,j) = s;
            X(j,i) = X(i,j);
        end;
    end;
    e1 = max(max(abs((X-M)/M)));
    disp(['Dimension = ',num2str(m)]);
    disp(['    Max error in UD decomp in-place = ',num2str(e1)]);
    D = zeros(m);
    U = UD;
    for i=1:m;
        D(i,i) = UD(i,i);
        U(i,i) = 1;
    end;
    UD = UDinv(UD);
    %
    Di = zeros(m);
    Ui = UD;
    for i=1:m;
        Di(i,i) = UD(i,i);
        Ui(i,i) = 1;
    end;
    e2 = max(max(abs(D*Di-eye(m))));
    disp(['    Max err in D inversion in-place = ',num2str(e2)]);
    e3 = max(max(abs(U*Ui-eye(m))));
    disp(['    Max err in U inversion in-place = ',num2str(e3)]);
    M0 = M;
    M = SPDinvIP(M);
    e4 = max(max(abs(M0*M-eye(m))));
    disp(['    Max err in M inversion in-place = ',num2str(e4)]);
end;