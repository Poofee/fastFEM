% 2017-12-12
% by Poofee
% 验证新的group传输线法
% 这次采用的是双传输线法
% 寻找最优猜测值
% 这种方法将一个小电路进行隔离cir2.jpg

close all;
clear all;
I = 10;
R1 = 10;Y1 = 1/R1;
R2 = 20;Y2 = 1/R2;
R3 = 30;Y3 = 1/R3;
R4 = 40;Y4 = 1/R4;
R5 = 50;Y5 = 1/R5;
R6 = 60;Y6 = 1/R6;
Y = [Y1+Y2+Y3, -Y2, -Y3;
    -Y2,Y2+Y4+Y5,-Y4;
    -Y3,-Y4,Y3+Y4+Y6];
I1 = [I;0;0];
V = Y\I1;
Z1 = 50;Yz1 = 1/Z1;%39.32
Z2 = 50;Yz2 = 1/Z2;
Z3 = 50;Yz3 = 1/Z3;
U1 = [0;0;0];
U2 = [0;0;0];
Vi = [0;0;0];
Vr = [0;0;0];
Yz = [Yz1;Yz2;Yz3];

for i=1:20
% 入射到小三角
Yi = [Y2+Y3+Yz1+Yz3,-Y2-Yz1;
    -Y2-Yz1,Y2+Y4+Yz1+Yz2];
%     -Y3-Yz3,-Y4-Yz2,Y3+Y4+Yz2+Yz3];
Ii = 2*[Vi(1)*Yz1-Vi(3)*Yz3;
    Vi(2)*Yz2-Vi(1)*Yz1;];
%     Vi(3)*Yz3-Vi(2)*Yz2];
U2 = Yi\Ii;
Vr = [U2(1)-U2(2);
    U2(2)-0;
    0-U2(1)]-Vi;
% 反射到另一部分
Yr = [Y1+Y5+Yz1+Yz3,-Yz1,-Yz3;
    -Yz1,Y5+Yz1+Yz2,-Yz2;
    -Yz3,-Yz2,Y6+Yz2+Yz3];
Ir = [I;0;0]+2*[Vr(1)*Yz1-Vr(3)*Yz3;
    Vr(2)*Yz2-Vr(1)*Yz1;
    Vr(3)*Yz3-Vr(2)*Yz2];
U1 = Yr\Ir;
Vi =  [U1(1)-U1(2);
    U1(2)-U1(3);
    U1(3)-U1(1)] - Vr;
plot(i,U1(1),'ro');hold on;
plot(i,U1(2),'bo');hold on;
plot(i,U1(3),'ko');hold on;
plot(i,V(1),'r*');hold on;
plot(i,V(2),'b*');hold on;
plot(i,V(3),'k*');hold on;
er = [V-U1];
% if abs(mean(er)) < 1e-3
%     break;
% end
end
