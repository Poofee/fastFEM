% 2017-11-24
% by Poofee
% 主要测试MTM方法
%参考文献
% 《transmission line inspires a new distributed 
% algorithm to solve the nonlinear dynamical
% system of physical circuit》
% 曲线表明这两种方法是等价的
close all
u = 10;%电压源
r1 = 5;%电阻1
r2 = 10;%电阻2
u1 = 0;
i1 = 0;
u2 = 0;
i2 = 0;
z = 200;%可调节参数，传输线阻抗

for i=1:100
    %left
    u1 = (u/r1 + u2/z - i2)/(1/r1 + 1/z);
    i1 = u2/z - i2 - u1/z;
    plot(i,u1,'r*');
    hold on
    %right
    u2 = (u1/z - i1)/(1/r2 + 1/z);
    i2 = u1/z - i1 - u2/z;
   plot(i,u2,'bo');
end

vi = 0;
vr = 0;
for i=1:100
    u1 = (u/r1 + 2*vi/z)/(1/r1 + 1/z);
    vr = u1 - vi;
    plot(i,u1,'bo')
    u2 = (2*vr/z)/(1/r2 + 1/z);
    vi = u2 - vr;
    plot(i,u2,'r*')
end