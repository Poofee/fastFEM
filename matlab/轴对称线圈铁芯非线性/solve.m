% 2017-10-30
% by Poofee
% 本例程使用三角形分网求解一个轴对称电磁结构
% 实现方法参考颜威利《电气工程电磁场数值分析》P53-P60
% 采用牛顿迭代法进行非线性迭代
clear all
close all
fname = ['mesh.mphtxt'];

[X,Y,NL,Domain]=readcomsoltri(fname);

% % 绘制所有的单元
%  for i=1:length(NL)
%      plot(X(NL(i,[1 2 3 1])),Y(NL(i,[1 2 3 1])));
%      hold on;
%  end
%   axis equal

num_nodes = length(X);
num_elements = length(NL);

CE = zeros(3,3);      % CE --- 用来存储每个单元的系数矩阵
D = zeros(3,3);
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

%该模型里只有5个域，域1，2为空气，域5为线圈，域3，4为铁芯。
%这是一个线性问题。
J = zeros(num_elements,1);%计算电流密度矩阵
coildomain = find(Domain == 5);%寻找线圈区域的单元
J(coildomain) = 8e5;%设置线圈区域的电流密度，其他其余为0

mu0 = 4*pi*1e-7;%空气磁导率
mu = mu0*ones(num_elements,1);%保存每一个单元的磁导率，初始化为空气磁导率

% 查找非边界上的节点
freenodes = find(abs(X)>1e-8 & sqrt(X.^2+Y.^2)<0.05-1e-3);
% plot(X,Y,'r*');
% hold on
% plot(X(freenodes),Y(freenodes),'bo');
% axis equal
% 计算教材上P59页的y
ydot = zeros(num_elements,1);
for i=1:num_elements
    %两个点在坐标轴上,注意P59页公式是错误的，应当为x
    if XL(i,1)+XL(i,2)<1e-10 || XL(i,2)+XL(i,3)<1e-10 || XL(i,1)+XL(i,3)<1e-10
        ydot(i) = mean(XL(i,:));
    else
        ydot(i) = 1.5/(1/(XL(i,1)+XL(i,2))+1/(XL(i,1)+XL(i,3))+1/(XL(i,2)+XL(i,3)));
    end
end
% 非线性迭代
A = zeros(num_nodes,1);%每个节点的磁势，轴对称情形下
B = zeros(num_elements,1);
Domain3 = find(Domain==3);
Domain4 = find(Domain==4);
Domain34 = [Domain3;Domain4];
steps = 20;
tol = 1e-10;%收敛误差
for count = 1:steps
%     装配
    for i=1:num_elements
        if Domain(i) == 3 || Domain(i) == 4
            dvdB = getdvdB(B(i));
        else%线性区域
            dvdB = 0;
        end
        
        for row = 1:3
            for col = 1:3
                D(row,col) = dvdB/ydot(i)/ydot(i)/ydot(i)/AREA(i);
                D(row,col) = D(row,col)* sum((R(i,row)*R(i,:)+Q(i,row)*Q(i,:)).*A(NL(i,:))')/4/AREA(i);
                D(row,col) = D(row,col)* sum((R(i,col)*R(i,:)+Q(i,col)*Q(i,:)).*A(NL(i,:))')/4/AREA(i);
                CE(row,col) = D(row,col) + (R(i,row)*R(i,col)+Q(i,row)*Q(i,col))/4/AREA(i)/mu(i)/ydot(i);
                
                S(NL(i,row),NL(i,col)) = S(NL(i,row),NL(i,col)) + CE(row,col);
                F1(NL(i,row)) = F1(NL(i,row)) + D(row,col)*A(NL(i,col));
            end            
            F1(NL(i,row)) = F1(NL(i,row)) + J(i)*AREA(i)/3;            
        end        
    end
    A_old = A;
    A(freenodes) = S(freenodes,freenodes)\F1(freenodes);
    % 判断误差
    error = norm((A_old - A))/norm(A);
    if error < tol
        break;
    end
%     下一步迭代初始化  
    S = S - S;
    F1 = F1 - F1;
    %更新B
    Bx = sum(R.*A(NL),2)./AREA./ydot/2;
    By = sum(Q.*A(NL),2)./AREA./ydot/2;
    B = sqrt(Bx.*Bx + By.*By);
    %更新mu，只更新非线性区域
    mu(Domain34) = B(Domain34)./arrayfun(@getH,B(Domain34));
end

A(freenodes) = A(freenodes) ./ X(freenodes);


F = scatteredInterpolant(X,Y,A);
% figure
% F = scatteredInterpolant(X,Y,A(1:num_nodes));
tx = 0:1e-4:0.025;
ty = -0.012:1e-4:0.012;
[qx,qy] = meshgrid(tx,ty);
qz = F(qx,qy);
% mesh(qx,qy,qz);
% surf(qx,qy,qz);
figure
subplot(1,2,2);hold on
title('MATLAB');
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
tx = 0:1e-4:0.025;
ty = -0.012:1e-4:0.012;
[qx,qy] = meshgrid(tx,ty);
qz = F(qx,qy);
% mesh(qx,qy,qz);
% surf(qx,qy,qz);
subplot(1,2,1);
contourf(qx,qy,qz,20);colorbar
title('COMSOL');
hold on
axis equal

% for i=1:length(Domain34)
%     plot(X(NL(Domain34,:)),Y(NL(Domain34,:)),'*')
%     hold on
% end
% axis equal