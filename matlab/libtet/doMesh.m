function [mesh] = doMesh(geofile)
% 20200419 by Poofee
% 调用gmsh对geo文件进行分网
tic
% 分网
fprintf('开始分网......\n');
cmd = ['gmsh.exe  -3 -format msh2 ',geofile];
[status,cmdout] = system(cmd);
mesh = load_gmsh2('model3d.msh');
fprintf('分网结束. 共 %d 个单元. 包括 %d 个节点, %d 个棱, %d 个三角形, %d 个四面体\n',...
    mesh.nbElm,mesh.nbNod,mesh.nbLines,mesh.nbTriangles,mesh.nbTets);
toc

% 绘制
% DisplayNodes(mesh.POS);
% 
% DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);
% 
% DisplayEdgeSB(100,mesh.TETS(:,1:4),mesh.POS);
end