% 2017-12-25
% by Poofee
% 验证单传输线法
% 寻找最优猜测值
close all
I = 10;%A
R1 = 10;
R2 = 40;
R3 = 30;
R4 = 40;
Y = [1/R1+1/R2+1/R4   -1/R2-1/R4;
   -1/R2-1/R4   1/R3+1/R2+1/R4];
I = [I;0];
V = Y\I;
Z1 = 40;
Vi = [0;0];
Vr = [0;0];
for i=1:50
    %incidence
    Y = [1/R1+1/R4+1/Z1   -1/R4;
   -1/R4   1/R3+1/R4 ];
   I1 = I + 2*Vi.*[1/Z1;0];
   V1 = Y\I1;
    %reflect
    Vr(1) = V1(1) - Vi(1);
    Y2 = [1/R2+1/Z1 -1/R2;
        -1/R2  1/R2+ 1/R3+1/(R1+R4)];
    I2 = 2*Vr.*[1/Z1;1/Z2];
    V2 = Y2\I2;
    Vi = V2 - Vr;
    plot(i,Vi(1),'ro');hold on;
    plot(i,Vi(2),'bo');hold on;
    plot(i,Vr(1),'r*');hold on;
    plot(i,Vr(2),'b*');hold on;
end
V1