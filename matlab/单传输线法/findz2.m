function step=findz2(z1,z2)
% 2017-12-05
% by Poofee
% 验证单传输线法
% 寻找最优猜测值
I = 10;%A
R1 = 1;
R2 = 2;
R3 = 3;
R4 = 4;
Y = [1/R1+1/R2+1/R4   -1/R2-1/R4;
   -1/R2-1/R4   1/R3+1/R2+1/R4];
I = [I;0];
V = Y\I;
Z1 = z1;
Z2 = z2;
Vi = [0;0];
Vr = [0;0];
for i=1:5000
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
    if(mean(abs(V-V1))<1e-6)
        step = i;
        break;
    end
%     plot(i,V1(1),'ro');hold on;
%     plot(i,V1(2),'bo');hold on;
%     plot(i,V2(1),'r*');hold on;
%     plot(i,V2(2),'b*');hold on;
end
step = i;
end