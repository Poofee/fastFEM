% 2017-12-05
% by Poofee
% 验证单传输线法
% 寻找最优猜测值
close all
I = 10;%A
R1 = 10;
R2 = 10;
R3 = 10;
R4 = 10;
Y = [1/R1+1/R2+1/R4   -1/R2-1/R4;
   -1/R2-1/R4   1/R3+1/R2+1/R4];
I = [I;0];
V = Y\I;
Z1 = 15;
Z2 = 15;
Vi = [0;0];
Vr = [0;0];
for i=1:50
    %incidence
    Y = [1/R1+1/R4+1/Z1   -1/R4;
   -1/R4   1/R3+1/R4+1/Z2];
   I1 = I + 2*Vi.*[1/Z1;1/Z2];
   V1 = Y\I1;
    %reflect
    Vr = V1 - Vi;
    Y2 = [1/R2+1/Z1 -1/R2;
        -1/R2  1/R2+1/Z2];
    I2 = 2*Vr.*[1/Z1;1/Z2];
    V2 = Y2\I2;
    Vi = V2 - Vr;
    plot(i,V1(1),'ro');hold on;
    plot(i,V1(2),'bo');hold on;
    plot(i,V2(1),'r*');hold on;
    plot(i,V2(2),'b*');hold on;
end
V1