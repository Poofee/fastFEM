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
cmd = 'gmsh.exe  -3 -format msh2 coil.geo';
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
NT = size(mesh.TETS,1);% ��Ԫ��Ŀ
Nn = size(mesh.POS,1);% �ڵ���Ŀ
% ���ɵ�Ԫ�����е��⣬˳���Ű���...
totalEdge = int32([mesh.TETS(:,[1 2]); mesh.TETS(:,[1 3]); mesh.TETS(:,[1 4]); ...
                   mesh.TETS(:,[2 3]); mesh.TETS(:,[2 4]); mesh.TETS(:,[3 4])]);
% �ԱߵĽڵ��Ž��д����ӵ͵���               
sortedTotalEdge = sort(totalEdge,2);
% edge��ȫ�ֵ��⣬��СNTx2�������ʼ��ͽ�����
% j������sortedTotalEdge�е�������edge�еı��
[edge,i2,j] = unique(sortedTotalEdge,'rows','legacy'); 
Ne = size(edge,1);% �����Ŀ
% elem2edge��ÿ����Ԫ��6�����ȫ�ֱ�ţ�˳��Ϊ[1,2][1,3][1,4][2,3][2,4][3,4]
elem2edge = uint32(reshape(j,NT,6));
direction = ones(6*NT,1);
idx = (totalEdge(:,1)>totalEdge(:,2));
direction(idx) = -1;
% elem2edgeSign�����Ǳ�Ŵӵ͵��ߵ���
elem2edgeSign = reshape(direction,NT,6);

%% ���㵥Ԫ����װ��
volume = zeros(mesh.nbTets,1);
base = mesh.nbPoints+mesh.nbLines+mesh.nbTriangles;
% ֻ���������ϰ벿�֣�����һ��6*6�ľ�����21��
rowAIndex = zeros(21*NT,1);% ������λ��
colAIndex = zeros(21*NT,1);% ������λ��
AA = zeros(21*NT,1);% �����Ӧ��Aֵ
% ����M����һ���Գƾ�������ȫ������
rowMIndex = zeros(24*NT,1);% ������λ��
colMIndex = zeros(24*NT,1);% ������λ��
MM = zeros(24*NT,1);% �����Ӧ��Mֵ

% �����Ҳ�����
% ��˹����
[lambda,w] = quadpts3(2);
nQuad = size(lambda,1);
locEdge = [1 2;1 3;1 4;2 3;2 4;3 4;];
bt = zeros(NT,6);

tmpIndex = 0;
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
    % �ⵥԪ���󣬴�СΪ6*6
    Ae = CurlW*CurlW'*volume(i)/mu0;
    % LM���Ӿ��󣬴�СΪ6*4
    Ge = [gradN(2,:)-gradN(1,:);
          gradN(3,:)-gradN(1,:);
          gradN(4,:)-gradN(1,:);
          gradN(3,:)-gradN(2,:);
          gradN(4,:)-gradN(2,:);
          gradN(4,:)-gradN(3,:);]*...
         [gradN(1,:);gradN(2,:);gradN(3,:);gradN(4,:)]'*volume(i)/20;
    % ����A��λ�þ���
    for cols=1:6
       for rows=cols:6
           tmpIndex = tmpIndex + 1;
           rowAIndex(tmpIndex) = elem2edge(i,rows);
           colAIndex(tmpIndex) = elem2edge(i,cols);
           AA(tmpIndex) = Ae(rows,cols);
       end
    end
    % ����M��λ�þ���
    for rows=1:6
       for cols=1:4
           Mindex = (i-1)*24+4*(rows-1)+cols;
           rowMIndex(Mindex) = elem2edge(i,rows);
           colMIndex(Mindex) = mesh.TETS(i,cols);
           MM(Mindex) = Ge(rows,cols);
       end
    end
    % �Ҳ�����
    for p = 1:nQuad
        % �����������
        pxyz = lambda(p,1)*mesh.POS(mesh.TETS(i,1),:)...
            +  lambda(p,2)*mesh.POS(mesh.TETS(i,2),:)...
            +  lambda(p,3)*mesh.POS(mesh.TETS(i,3),:)...
            +  lambda(p,4)*mesh.POS(mesh.TETS(i,4),:);
        % �����pxyz���ĵ����ܶ�
        Jp = [0,0,0];        
        if mesh.ELE_TAGS(base+i,2) == CoilTag
            % ����Ƕ�
            costheta = pxyz(1)/norm([pxyz(1),pxyz(2)]);
            if pxyz(2) > 0
                sintheta = sqrt(1-costheta*costheta);
            else
                sintheta = -sqrt(1-costheta*costheta);
            end
            Jp = Js*[-sintheta,costheta,0];
        end
        for k=1:6
            a1 = locEdge(k,1);
            a2 = locEdge(k,2);
            phi_k = lambda(p,a1)*gradN(a2,:)-lambda(p,a2)*gradN(a1,:);
            rhs = dot(phi_k,Jp);
            bt(i,k) = bt(i,k) + w(p)*rhs;
        end
    end
end

% ����ϡ�����
% �Խ����ϵ�Ԫ����������Ԫ��
diagAIndex = (rowAIndex == colAIndex);
upperAIndex = ~diagAIndex;
% ϡ�����A�Ĵ�СΪNe*Ne
A = sparse(rowAIndex(diagAIndex),colAIndex(diagAIndex),AA(diagAIndex),Ne,Ne);
AUpper = sparse(rowAIndex(upperAIndex),colAIndex(upperAIndex),AA(upperAIndex),Ne,Ne);
A = A + AUpper + AUpper';% ���ݶԳ��Բ�ȫ�����
M = sparse(rowMIndex,colMIndex,MM,Ne,Nn);
% ����f
bt = bt.*repmat(volume,1,6);
f = accumarray(elem2edge(:),bt(:),[Ne 1]);

%% ʩ�ӱ߽�����
% һ������ʱ�򣬲���Ҫ�ر�����ñ߽���������Ϊ��һ��Ĭ�ϵļ��α߽磬��ôΪ��
% ���ñ߽磬����Ҫ�ҳ�����߽硣��ǰ��˼·���ü��γߴ����жϣ�ȱ����������֪��
% �߽����״�ͳߴ磬�ǳ������㡣���ǣ�����һ�£����һ����Ԫ�Ǳ߽磬��ά��ʱ��
% �����߽�ֻ������һ��������Ԫ�У���ô�Ϳ��Բ��ҳ���һ�εıߣ��Ǿ��Ǳ߽��ˡ�
% ��ά��˼·Ҳ��һ���ģ�ֻ����Ҫ���ҵ����浥Ԫ��
% �����������������
allFace = [mesh.TETS(:,[2 4 3]);mesh.TETS(:,[1 3 4]);mesh.TETS(:,[1 4 2]);mesh.TETS(:,[1 2 3])];

Nfall = length(allFace);% �����Ŀ
% ���ظ���face
[face, i2, j] = unique(sort(allFace,2),'rows','legacy');

i1(j(Nfall:-1:1)) = Nfall:-1:1;% ֻ�г���һ�εı�ŲŲ���� 
i1 = i1';
bdFlag = zeros(Nfall,1,'uint8');
bdFaceidx = i1(i1==i2);% �߽��ϵ�face������
bdFlag(bdFaceidx) = 1;% Dirichlet�߽�����
bdFlag = reshape(bdFlag,NT,4);% ת��Ϊ��Ӧ�ĵ�Ԫ��ʽ
% ���ұ߽���ͽڵ�
isBdEdge = false(Ne,1);
isBdNode = false(Nn,1);
isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
bdEdge = edge(isBdEdge,:);% ����ά����ת��Ϊ����
isBdNode(bdEdge) = true;% ��Ǳ߽�ڵ�
% δ֪�ı���
freeEdge = find(~isBdEdge);
freeNode = find(~isBdNode);
% ����Neumann�߽�����

% ����Dirichlet�߽�����
u = zeros(Ne,1);
if ~isempty(bdEdge)
    u(isBdEdge) = 0;% �̶��߽�����
%     f = f - (A - M)*u;
    f(isBdEdge) = u(isBdEdge);
end
g0 = -M'*u;
%% ���
% reduce to free dofs
A  = A(freeEdge,freeEdge);
M  = M(freeEdge,freeNode);
f = f(freeEdge);
g0 = g0(freeNode);
Ni = size(freeNode,1);
bigA = [A M;...
        M' sparse(Ni,Ni)];
tmp = bigA\[f;g0];

%% ��������