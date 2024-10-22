% 2018-03-04
% by Poofee
% 本程序是为了验证coupled transmission line的时候，
% 传输线导纳是否可以设置为一个对称非对角矩阵。
% 这个电路是三端
% 电路图采用cir3.jpg
% 传输线导纳当中含地
close all
clear all
I1 = 10;%A
R1 = 10;
R2 = 20;
R3 = 30;
R4 = 40;
R5 = 50;
R6 = -10;
R7 = 10;

% 求解电路解
Y = [1/R1+1/R2, -1/R2,  0,  0;
        -1/R2,  1/R2+1/R3+1/R4,-1/R3,-1/R4;
        0,  -1/R3,  1/R3+1/R5+1/R6, -1/R6;
        0,  -1/R4,  -1/R6,  1/R4+1/R6+1/R7;];
I = [I1;0;0;0];
V = Y\I;

% 迭代需要的变量
U1 = [10;10;10;10];
i1 = [0;0;0];
U2 = [10;10;0];
i2 = [0;0;0];
a = [0;0;0];
b = [0;0;0];
% 左边的导纳矩阵
Y = [1/R1+1/R2,    -1/R2,  0,           0;
        -1/R2,            1/R2,    0,           0;
        0,                   0,          1/R5,     0;
        0,                   0,          0,           1/R7;];
% 设置传输线导纳，设为真值
R3c = 30;
R4c = 40;
R6c = 50;
YTL = [1/R3c+1/R4c+1/10000,        -1/R3c,             -1/R4c;
          -1/R3c,                1/R3c+1/R6c+0,     -1/R6c;
          -1/R4c,                -1/R6c,              1/R6c+1/R4c+0;]; 
YTL = 1* YTL;
Y(2:4,2:4) = Y(2:4,2:4) + YTL;
% 右边的导纳矩阵
Y2 = [1/R3+1/R4,        -1/R3,             -1/R4;
          -1/R3,                1/R3+1/R6,     -1/R6;
          -1/R4,                -1/R6,              1/R6+1/R4;];   
Y2 = Y2 + YTL;
for i=1:5000
%     左边的电流源电流
%     a(1) = U2(1)/Z1-i2(1);
%     a(2) = U2(2)/Z2-i2(2);
    aa = a;
    bb = b;
    a = 2*YTL*U2 - bb;
%     右边的电流源电流
%     b(1) = U1(1)/Z1-i1(1);
%     b(2) = U1(2)/Z2-i1(2);
    b = 2*YTL*U1(2:4) - aa;
    %incidence   
    IL = I;
    IL(2:4) = IL(2:4) + a;    
    %reflect
    I2 = b;
    U1 = Y\IL;
    U2 = Y2\I2;%右边的电路跟左边的地不一样，降阶矩阵    
    
%     i1(1) = a(1)-U1(1)/Z1;
%     i1(2) = a(2)-U1(2)/Z2;
%     i1 = a - YTL*U1(2:4);
    
%     i2(1) = b(1)-U2(1)/Z1;
%     i2(2) = b(2)-U2(2)/Z2;
%     i2 = b - YTL*U2;
    
    plot(i,U1(2),'ro');hold on;
    plot(i,U1(3),'bo');hold on;
    plot(i,U1(4),'go');hold on;
%     plot(i,U2(1),'r*');hold on;
%     plot(i,U2(2),'b*');hold on;
%     plot(i,U2(3),'g*');hold on;
    plot(i,V(2),'r.-');hold on;
    plot(i,V(3),'b.-');hold on;
    plot(i,V(4),'g.-');hold on;
    if (mean(abs(U1-V)) > 1e-8) && (i > 2)
%         i
        break;
    end
end
