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

% plot3(mesh.POS(:,1),mesh.POS(:,2),mesh.POS(:,3),'*','MarkerSize',24);
% for i=1:mesh.nbNod
%     text(mesh.POS(i,1),mesh.POS(i,2),mesh.POS(i,3),num2str(i),'Color',[1 0 0],'FontSize',20);
% end
% axis equal
% hold on
for i=1:mesh.nbTets
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(i,1:4),3);
    tmp = [ones(4,1),X,Y,Z];
    a = det(tmp)/6;
    if a < 0
       disp(['警告：单元',num2str(i),'编号不规范']) 
    end
    volume(i) = a;
%     if mesh.ELE_TAGS(base+i,2)==80
%         Lines = [1 2 3 4 4 1;...
%              2 3 4 1 2 3];
%         line(X(Lines),Y(Lines),Z(Lines),'Color',[0 0 0]);
%         fill3(X([1,2,3]),Y([1,2,3]),Z([1,2,3]),[0 0 1],'facealpha',0.5);hold on
%         fill3(X([1,2,4]),Y([1,2,4]),Z([1,2,4]),[0 1 1],'facealpha',0.5);
%         fill3(X([2,3,4]),Y([2,3,4]),Z([2,3,4]),[0 1 1],'facealpha',0.5);
%         fill3(X([1,3,4]),Y([1,3,4]),Z([1,3,4]),[0 1 1],'facealpha',0.5);
%         text(X,Y,Z,{'1','2','3','4'},'Color',[1 0 0],'FontSize',24);
%         axis equal
%         drawnow
%     end
end



