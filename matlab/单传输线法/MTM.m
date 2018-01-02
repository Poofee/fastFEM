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
Z1 = 10;
Z2 = 100;
U1 = [0;0];
i1 = [0;0];
U2 = [0;0];
i2 = [0;0];
for i=1:50
    %incidence
    Y = [1/R1+1/R4+1/Z1   -1/R4;
   -1/R4   1/R3+1/R4+1/Z2];
   I1 = I + [U2(1)/Z1-i2(1);U2(2)/Z2-i2(2)];
   U1 = Y\I1;
   i1(1) = U2(1)/Z1-i2(1)-U1(1)/Z1;
   i1(2) = U2(2)/Z2-i2(2)-U1(2)/Z2;
    %reflect
    Y2 = [1/R2+1/Z1 -1/R2;
        -1/R2  1/R2+1/Z2];
    I2 = [U1(1)/Z1-i1(1);U1(2)/Z2-i1(2)];
    U2 = Y2\I2;
    i2(1) = U1(1)/Z1-i1(1)-U2(1)/Z1;
    i2(2) = U1(2)/Z2-i1(2)-U2(2)/Z2;
    plot(i,U1(1),'ro');hold on;
    plot(i,U1(2),'bo');hold on;
    plot(i,U2(1),'r*');hold on;
    plot(i,U2(2),'b*');hold on;
end
U1-V