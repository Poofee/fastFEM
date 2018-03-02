% 2018-03-02
% by Poofee
% 本程序是为了验证coupled transmission line的时候，
% 传输线导纳是否可以设置为一个对称非对角矩阵。
% 电路图采用cir.png
close all
clear all
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
Z2 = 10;
U1 = [0;0];
i1 = [0;0];
U2 = [0;0];
i2 = [0;0];
a = [0;0];
b = [0;0];
% 左边的导纳矩阵
Y = [1/R1+1/R4   -1/R4;
        -1/R4   1/R3+1/R4];
 YL = [1/10+1/10,-1/10;
           -1/10,     1/10+1/10];
Y = Y + YL;
% 右边的导纳矩阵
Y2 = [1/R2 -1/R2;
        -1/R2  1/R2];   
YR = YL;
Y2 = Y2 + YR;
for i=1:50
%     左边的电流源电流
%     a(1) = U2(1)/Z1-i2(1);
%     a(2) = U2(2)/Z2-i2(2);
    a = YL*U2 - i2;
%     右边的电流源电流
%     b(1) = U1(1)/Z1-i1(1);
%     b(2) = U1(2)/Z2-i1(2);
    b = YR*U1 - i1;
    %incidence   
    I1 = I + a;    
    %reflect
    I2 = b;
    U1 = Y\I1;
    U2 = Y2\I2;    
    
%     i1(1) = a(1)-U1(1)/Z1;
%     i1(2) = a(2)-U1(2)/Z2;
    i1 = a - YR*U1;
    
%     i2(1) = b(1)-U2(1)/Z1;
%     i2(2) = b(2)-U2(2)/Z2;
    i2 = b - YL*U2;
    
    plot(i,U1(1),'ro');hold on;
    plot(i,U1(2),'bo');hold on;
    plot(i,U2(1),'r*');hold on;
    plot(i,U2(2),'b*');hold on;
    if abs(U1-V) < 1e-3
        i
        break;
    end
end
