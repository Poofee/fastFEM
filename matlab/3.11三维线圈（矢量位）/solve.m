% 20191025 by Poofee
% ����ʸ���ⵥԪ������ά���ų��������
% ���ģ��Ϊ��Ȧ�����Ĵų�

% ��Ȧ��һЩ������
% ��Ȧ��ѹ��25 V
% ��Ȧ������1000
% ��Ȧ�ھ���16 mm
% ��Ȧ�⾶��24 mm
% ��Ȧ�߶ȣ�46 mm
close all
clear all
tic
%% ����gmsh����
fprintf('��ʼ����......\n');
cmd = ['gmsh.exe  -3 -format msh2 coil.geo'];
[status,cmdout] = system(cmd);
mesh = load_gmsh2('coil.msh');
fprintf('��������. �� %d ����Ԫ. ���� %d ���ڵ�, %d ����, %d ��������, %d ��������\n',...
    mesh.nbElm,mesh.nbPoints,mesh.nbLines,mesh.nbTriangles,mesh.nbTets);
toc

%% ��ȡCOMSOL����
% mesh = ReadCOMSOL3D('coil.mphtxt');
%% ����
% DisplayNodes(mesh.POS);

% DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);

% DisplayEdgeSB(100,mesh.TETS(:,1:4),mesh.POS);

%% �������
Ucoil = 25;%��Ȧ��ѹ
Ncoil = 1000;%��Ȧ����
Rc = 10;%��Ȧ����
Scoil = 0.46 * (0.24 - 0.16);%��Ȧ���
Js = Ncoil*Ucoil/Rc/Scoil;%��Ȧ�����ܶ�
mu0 = 4*pi*1e-7;%�����ŵ���
AirTag = 80;
CoilTag = 1;

%% ���㵥Ԫ����
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
       disp(['���棺��Ԫ',num2str(i),'��Ų��淶']) 
    end
end

%% �����Ҳ�����

%% ʩ�ӱ߽�����

%% ���


