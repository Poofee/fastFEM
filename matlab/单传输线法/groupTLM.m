% 2017-12-07
% by Poofee
% 验证group传输线法
% 寻找最优猜测值
% 这种方法将一个小电路进行隔离
close all;
clear all;
I = 10;
R1 = 10;Y1 = 1/R1;
R2 = 10;Y2 = 1/R2;
R3 = 10;Y3 = 1/R3;
R4 = 10;Y4 = 1/R4;
R5 = 10;Y5 = 1/R5;
R6 = 10;Y6 = 1/R6;
Y = [Y1+Y3+Y4, -Y3, -Y4;
    -Y3,Y2+Y3+Y5+Y6,-Y5-Y6;
    -Y4,-Y5-Y6,Y4+Y5+Y6];
I1 = [I;0;0];
V = Y\I1;
Z1 = 10;Yz1 = 1/Z1;
Z2 = 10;Yz2 = 1/Z2;
Z3 = 10;Yz3 = 1/Z3;
U1 = [0;0;0];
U2 = [0;0;0];
Vi = [0;0;0];
Vr = [0;0;0];
Yz = [Yz1;Yz2;Yz3];

for i=1:100
% 入射到小三角
Yi = [Y3+Y4+Yz1,-Y3,-Y4;
    -Y3,Y3+Y5+Yz2,-Y5;
    -Y4,-Y5,Y4+Y5+Yz3];
Ii = 2*Vi.*Yz;
U2 = Yi\Ii;
Vr = U2 - Vi;
% 反射到另一部分
Yr = [Y1+Yz1,0,0;
    0,Y2+Yz2+Y6,-Y6;
    0,-Y6,Y6+Yz3];
Ir = [I;0;0]+2*Vr.*Yz;
U1 = Yr\Ir;
Vi = U1 - Vr;
plot(i,U1(1),'ro');hold on;
plot(i,U1(2),'bo');hold on;
plot(i,U1(3),'yo');hold on;
plot(i,U2(1),'r*');hold on;
plot(i,U2(2),'b*');hold on;
plot(i,U2(3),'y*');hold on;
end
[V-U1,V-U2]