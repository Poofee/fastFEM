% 20190726 By Poofee
% ��ά���ų����
% ģ��ΪFLUX����


tic
% ����
fprintf('��ʼ����......\n');
cmd = ['gmsh.exe  -3 -format msh2 model3d.geo'];
[status,cmdout] = system(cmd);
mesh = load_gmsh2('model3d.msh');
fprintf('��������. �� %d ����Ԫ. ���� %d ���ڵ�, %d ����, %d ��������, %d ��������\n',...
    mesh.nbElm,mesh.nbNod,mesh.nbLines,mesh.nbTriangles,mesh.nbTets);
toc

% ����
% DisplayNodes(mesh.POS);

DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);