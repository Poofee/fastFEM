% 2017-10-29
% by Poofee
% 本例程使用三角形分网求解线圈产生的静磁场
% 实现方法参考颜威利《电气工程电磁场数值分析》P46-P50
clear all
close all
fname = ['coilmesh.mphtxt'];

[X,Y,NL,Domain]=readcomsoltri(fname);

% % 绘制所有的单元
%  for i=1:length(NL)
%      plot(X(NL(i,[1 2 3 1])),Y(NL(i,[1 2 3 1])));
%      hold on;
%  end
%   axis equal

num_nodes = length(X);
num_elements= length(NL);

CE = zeros(3,3);      % CE --- 用来存储每个单元的系数矩阵
S = zeros(num_nodes,num_nodes);%全局矩阵
F1 = zeros(num_nodes,1);

%AREA = zeros(num_elements,1);%三角形面积
P = zeros(num_elements,3);
Q = zeros(num_elements,3);
R = zeros(num_elements,3);

XL = X(NL);
YL = Y(NL);
Q(:,1) = YL(:,2) - YL(:,3);
Q(:,2) = YL(:,3) - YL(:,1);
Q(:,3) = YL(:,1) - YL(:,2);

R(:,1) = XL(:,3) - XL(:,2);
R(:,2) = XL(:,1) - XL(:,3);
R(:,3) = XL(:,2) - XL(:,1);

P(:,1) = XL(:,2).*YL(:,3) - XL(:,3).*YL(:,2);
P(:,2) = XL(:,3).*YL(:,1) - XL(:,1).*YL(:,3);
P(:,3) = XL(:,1).*YL(:,2) - XL(:,2).*YL(:,1);

AREA = 0.5 * (Q(:,1).*R(:,2) - Q(:,2).*R(:,1));%三角形面积

%该模型里只有两个域，域1为空气，域2为线圈。
%这是一个线性问题。
J = zeros(num_elements,1);%计算电流密度矩阵
coildomain = find(Domain == 2);%寻找线圈区域的单元
J(coildomain) = 8e5;%设置线圈区域的电流密度，其他其余为0

mu0 = 4*pi*1e-7;%空气磁导率

for i=1:num_elements
    for row = 1:3
        for col = 1:3
            CE(row,col) = (R(i,row)*R(i,col)+Q(i,row)*Q(i,col))/4/AREA(i)/mu0;
            S(NL(i,row),NL(i,col)) = S(NL(i,row),NL(i,col)) + CE(row,col);
        end
        F1(NL(i,row)) = F1(NL(i,row)) + J(i)*AREA(i)/3;
    end
end
%查找边界点
A = zeros(num_nodes,1);
freenodes = find(abs(X)~=1 & abs(Y) ~=1);
A(freenodes) = S(freenodes,freenodes)\F1(freenodes);

FF = scatteredInterpolant(X,Y,A);
% figure
% F = scatteredInterpolant(X,Y,A(1:num_nodes));
tx = -1:1e-2:1;
ty = -1:1e-2:1;
[qx,qy] = meshgrid(tx,ty);
qz = FF(qx,qy);
% mesh(qx,qy,qz);
% surf(qx,qy,qz);
figure
subplot(1,2,2);hold on
title('MATLAB','FontName','Times New Roman','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15);
contourf(qx,qy,qz,20);colorbar
axis equal

% 读入COMSOL的数据文件，进行云图绘制
fp = fopen('comsoldata.txt','r');
% 读取前面没用的9行
for i = 1:9
    fgets(fp);
end
% 读取三列数据
comsoldata = fscanf(fp,'%lf %lf %lf\n',[ 3 num_nodes]);
comsoldata = comsoldata';
fclose(fp);

F = scatteredInterpolant(comsoldata(:,1),comsoldata(:,2),comsoldata(:,3));
% figure
% F = scatteredInterpolant(X,Y,A(1:num_nodes));
tx = -1:1e-2:1;
ty = -1:1e-2:1;
[qx,qy] = meshgrid(tx,ty);
qz = F(qx,qy);
% mesh(qx,qy,qz);
% surf(qx,qy,qz);
subplot(1,2,1);
contourf(qx,qy,qz,20);colorbar
title('COMSOL','FontName','Times New Roman','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15);
hold on
axis equal