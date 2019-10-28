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

%% ��������Ŀ
NT = size(mesh.TETS,1);
totalEdge = int32([mesh.TETS(:,[1 2]); mesh.TETS(:,[1 3]); mesh.TETS(:,[1 4]); ...
                   mesh.TETS(:,[2 3]); mesh.TETS(:,[2 4]); mesh.TETS(:,[3 4])]);
% �Ա߽�������               
sortedTotalEdge = sort(totalEdge,2);
% edge��ȫ�ֵ��⣬��СNTx2�������ʼ��ͽ�����
[edge,i2,j] = unique(sortedTotalEdge,'rows','legacy'); 
elem2edge = uint32(reshape(j,NT,6));% ÿ����Ԫ��6�����ȫ�ֱ��
direction = ones(6*NT,1);
idx = (totalEdge(:,1)>totalEdge(:,2));
direction(idx) = -1;
elem2edgeSign = reshape(direction,NT,6);

%% ���㵥Ԫ����
volume = zeros(mesh.nbTets,1);
base = mesh.nbPoints+mesh.nbLines+mesh.nbTriangles;

for i=1:mesh.nbTets
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(i,1:4),3);
    tmp = [ones(4,1),X,Y,Z];
    % ���㵥Ԫ���
    volume(i) = det(tmp)/6;
    if volume(i) < 0
       volume(i) = - volume(i);
       disp(['���棺��Ԫ',num2str(i),'��Ų��淶']) 
    end
    gradN = zeros(4,3);
    gradN(1,:) = dTetraNodalBasis(1,X,Y,Z,0,0,0);
    gradN(2,:) = dTetraNodalBasis(2,X,Y,Z,0,0,0);
    gradN(3,:) = dTetraNodalBasis(3,X,Y,Z,0,0,0);
    gradN(4,:) = dTetraNodalBasis(4,X,Y,Z,0,0,0);
    % ���� \nabla \times \mathbf{w}_{i}
    CurlW = RotWBasis(X,Y,Z,0,0,0);
    % �ⵥԪ����
    Ae = CurlW*CurlW'*volume(i)/mu0;
    % LM���Ӿ���
    Ge = [gradN(2,:)-gradN(1,:);
          gradN(3,:)-gradN(1,:);
          gradN(4,:)-gradN(1,:);
          gradN(3,:)-gradN(2,:);
          gradN(4,:)-gradN(2,:);
          gradN(4,:)-gradN(3,:);]*...
          [gradN(1,:);gradN(2,:);gradN(3,:);gradN(4,:)]'*volume(i)/20;
end

%% �����Ҳ�����
% ��˹����
[lambda,w] = quadpts3(2);
nQuad = size(lambda,1);
for i=1:mesh.nbTets
    % ��������
    if mesh.ELE_TAGS(base+i,2) == CoilTag
        
    end
end
%% ʩ�ӱ߽�����

%% ���


