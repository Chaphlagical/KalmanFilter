%
% plots pYArctanaX
%
close all;
clear all;
n       = 0;
for log10a=-1:.02:1,
    n       = n + 1;
    a(n)    = 10^log10a;
    k       = 0;
    for y=-pi/2:pi/100:pi/2,
        k       = k + 1;
        if n==1
            Y(k)    = y;
        end;
        pY(n,k) = pYArctanaX(y,a(n));
    end;
end;
[X,Y]       = meshgrid(-10:.2:10);
[rows,cols] = size(X);
C           = zeros(rows,cols);
figure;
meshc(X,Y,pY,C);
axis([-10 10 -10 10 0 1.3]);
axis off;
view(60,40);
figure;
a1      = 1;
a3     = length(pY(:,1));
a2      = (a1+a3)/2;
plot(Y,pY(a1,:),'k-',Y,pY(a2,:),'k--',Y,pY(a3,:),'k-.','LineWidth',1.5);
xlabel('y','FontSize',14);
ylabel('p_y(y)','FontSize',14);
legend(['a=',num2str(a(a1))],['a=',num2str(a(a2))],['a=',num2str(a(a3))]);
title('PROBABILITY DENSITIES OF arctan( a X ), X ZERO-MEAN UNIT GAUSSIAN','FontSize',12);

        