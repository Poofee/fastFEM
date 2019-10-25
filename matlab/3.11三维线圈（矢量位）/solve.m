% 20191025 by Poofee
% 采用矢量棱单元法对三维静磁场进行求解
% 求解模型为线圈产生的磁场

% 线圈的一些参数：
% 线圈电压：25 V
% 线圈匝数：1000
% 线圈内径：16 mm
% 线圈外径：24 mm
% 线圈高度：46 mm

tic
% 分网
fprintf('开始分网......\n');
cmd = ['gmsh.exe  -3 -format msh2 coil.geo'];
[status,cmdout] = system(cmd);
mesh = load_gmsh2('coil.msh');
fprintf('分网结束. 共 %d 个单元. 包括 %d 个节点, %d 个棱, %d 个三角形, %d 个四面体\n',...
    mesh.nbElm,mesh.nbNod,mesh.nbLines,mesh.nbTriangles,mesh.nbTets);
toc

% 绘制
DisplayNodes(mesh.POS);

DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);

% DisplayEdgeSB(100,mesh.TETS(:,1:4),mesh.POS);

