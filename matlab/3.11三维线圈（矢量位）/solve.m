% 20191025 by Poofee
% ����ʸ���ⵥԪ������ά���ų��������
% ���ģ��Ϊ��Ȧ�����Ĵų�

% ��Ȧ��һЩ������
% ��Ȧ��ѹ��25 V
% ��Ȧ������1000
% ��Ȧ�ھ���16 mm
% ��Ȧ�⾶��24 mm
% ��Ȧ�߶ȣ�46 mm

tic
% ����
fprintf('��ʼ����......\n');
cmd = ['gmsh.exe  -3 -format msh2 coil.geo'];
[status,cmdout] = system(cmd);
mesh = load_gmsh2('coil.msh');
fprintf('��������. �� %d ����Ԫ. ���� %d ���ڵ�, %d ����, %d ��������, %d ��������\n',...
    mesh.nbElm,mesh.nbNod,mesh.nbLines,mesh.nbTriangles,mesh.nbTets);
toc

% ����
DisplayNodes(mesh.POS);

DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);

% DisplayEdgeSB(100,mesh.TETS(:,1:4),mesh.POS);

