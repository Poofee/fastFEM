% 2017-11-24
% by Poofee
% ��Ҫ����MTM����
%�ο�����
% ��transmission line inspires a new distributed 
% algorithm to solve the nonlinear dynamical
% system of physical circuit��
% ���߱��������ַ����ǵȼ۵�
close all
u = 10;%��ѹԴ
r1 = 5;%����1
r2 = 10;%����2
u1 = 0;
i1 = 0;
u2 = 0;
i2 = 0;
u1a = 0;
u2a = 0;
i1a = 0;
i2a = 0;
z = 11;% �ɵ��ڲ������������迹
t = zeros(1,20);
% for z = 5:5:20
    u1 = 0;
    i1 = 0;
    u2 = 0;
    i2 = 0;
    for i=1:50
        %left
        u1a = u1;
        i1a = i1;
        u1 = (u/r1 + u2a/z - i2a)/(1/r1 + 1/z);
        i1 = u2a/z - i2a - u1/z;
        plot(i,u1,'b*');
        t(i) = u1;
        hold on
        %right
        u2a = u2;
        i2a = i2;
        u2 = (u1a/z - i1a)/(1/r2 + 1/z);
        i2 = u1a/z - i1a - u2/z;
       plot(i,u2,'bo');
    end
    plot(t,'*-')
%     legend(['z=',num2str(z)])
% end

vi = 0;
vr = 0;
z=10
for i=1:50
    u1 = (u/r1 + 2*vi/z)/(1/r1 + 1/z);
    vr = u1 - vi;
    plot(i,u1,'ro')
    u2 = (2*vr/z)/(1/r2 + 1/z);
    vi = u2 - vr;
    plot(i,u2,'r*')
end