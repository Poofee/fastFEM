% 20191025 by Poofee
% 采用矢量棱单元法对三维静磁场进行求解
% 求解模型为线圈产生的磁场

% 线圈的一些参数：
% 线圈电压：25 V
% 线圈匝数：1000
% 线圈内径：16 mm
% 线圈外径：24 mm
% 线圈高度：46 mm
close all
clear all
tic
%% 调用gmsh分网
fprintf('开始分网......\n');
cmd = ['gmsh.exe  -3 -format msh2 coil.geo'];
[status,cmdout] = system(cmd);
mesh = load_gmsh2('coil.msh');
fprintf('分网结束. 共 %d 个单元. 包括 %d 个节点, %d 个棱, %d 个三角形, %d 个四面体\n',...
    mesh.nbElm,mesh.nbPoints,mesh.nbLines,mesh.nbTriangles,mesh.nbTets);
toc

%% 读取COMSOL分网
% mesh = ReadCOMSOL3D('coil.mphtxt');
%% 绘制
% DisplayNodes(mesh.POS);

% DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);

% DisplayEdgeSB(100,mesh.TETS(:,1:4),mesh.POS);

%% 物理参数
Ucoil = 25;%线圈电压
Ncoil = 1000;%线圈匝数
Rc = 10;%线圈电阻
Scoil = 0.46 * (0.24 - 0.16);%线圈面积
Js = Ncoil*Ucoil/Rc/Scoil;%线圈电流密度
mu0 = 4*pi*1e-7;%空气磁导率
AirTag = 80;
CoilTag = 1;

%% 计算单元矩阵
volume = zeros(mesh.nbTets,1);
base = mesh.nbPoints+mesh.nbLines+mesh.nbTriangles;

for i=1:mesh.nbTets
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(i,1:4),3);
    tmp = [ones(4,1),X,Y,Z];
    volume(i) = det(tmp)/6;
    if volume(i) < 0
       volume(i) = - volume(i);
       disp(['警告：单元',num2str(i),'编号不规范']) 
    end
end

%% 计算右侧向量

%% 施加边界条件

%% 求解


