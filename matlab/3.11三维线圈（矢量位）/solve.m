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
cmd = 'gmsh.exe  -3 -format msh2 coil.geo';
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

%% 计算棱数目
NT = size(mesh.TETS,1);% 单元数目
Nn = size(mesh.POS,1);% 节点数目
% 生成单元中所有的棱，顺序编号按照...
totalEdge = int32([mesh.TETS(:,[1 2]); mesh.TETS(:,[1 3]); mesh.TETS(:,[1 4]); ...
                   mesh.TETS(:,[2 3]); mesh.TETS(:,[2 4]); mesh.TETS(:,[3 4])]);
% 对边的节点编号进行处理，从低到高               
sortedTotalEdge = sort(totalEdge,2);
% edge，全局的棱，大小NTx2，棱的起始点和结束点
% j，保存sortedTotalEdge中的数据在edge中的编号
[edge,i2,j] = unique(sortedTotalEdge,'rows','legacy'); 
Ne = size(edge,1);% 棱的数目
% elem2edge，每个单元的6条棱的全局编号，顺序为[1,2][1,3][1,4][2,3][2,4][3,4]
elem2edge = uint32(reshape(j,NT,6));
direction = ones(6*NT,1);
idx = (totalEdge(:,1)>totalEdge(:,2));
direction(idx) = -1;
% elem2edgeSign，不是编号从低到高的棱
elem2edgeSign = reshape(direction,NT,6);

%% 计算单元矩阵及装配
volume = zeros(mesh.nbTets,1);
base = mesh.nbPoints+mesh.nbLines+mesh.nbTriangles;
% 只保存矩阵的上半部分，对于一个6*6的矩阵，有21个
rowAIndex = zeros(21*NT,1);% 保存行位置
colAIndex = zeros(21*NT,1);% 保存列位置
AA = zeros(21*NT,1);% 保存对应的A值
% 由于M不是一个对称矩阵，所以全部保存
rowMIndex = zeros(24*NT,1);% 保存行位置
colMIndex = zeros(24*NT,1);% 保存列位置
MM = zeros(24*NT,1);% 保存对应的M值

% 计算右侧向量
% 高斯积分
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
    % 计算单元体积
    volume(i) = det(tmp)/6;
    if volume(i) < 0
       volume(i) = - volume(i);
       disp(['警告：单元',num2str(i),'编号不规范']) 
    end
    gradN = zeros(4,3);
    gradN(1,:) = dTetraNodalBasis(1,X,Y,Z,0,0,0);
    gradN(2,:) = dTetraNodalBasis(2,X,Y,Z,0,0,0);
    gradN(3,:) = dTetraNodalBasis(3,X,Y,Z,0,0,0);
    gradN(4,:) = dTetraNodalBasis(4,X,Y,Z,0,0,0);
    % 计算 \nabla \times \mathbf{w}_{i}
    CurlW = RotWBasis(X,Y,Z,0,0,0);
    % 棱单元矩阵，大小为6*6
    Ae = CurlW*CurlW'*volume(i)/mu0;
    % LM算子矩阵，大小为6*4
    Ge = [gradN(2,:)-gradN(1,:);
          gradN(3,:)-gradN(1,:);
          gradN(4,:)-gradN(1,:);
          gradN(3,:)-gradN(2,:);
          gradN(4,:)-gradN(2,:);
          gradN(4,:)-gradN(3,:);]*...
         [gradN(1,:);gradN(2,:);gradN(3,:);gradN(4,:)]'*volume(i)/20;
    % 保存A的位置矩阵
    for cols=1:6
       for rows=cols:6
           tmpIndex = tmpIndex + 1;
           rowAIndex(tmpIndex) = elem2edge(i,rows);
           colAIndex(tmpIndex) = elem2edge(i,cols);
           AA(tmpIndex) = Ae(rows,cols);
       end
    end
    % 保存M的位置矩阵
    for rows=1:6
       for cols=1:4
           Mindex = (i-1)*24+4*(rows-1)+cols;
           rowMIndex(Mindex) = elem2edge(i,rows);
           colMIndex(Mindex) = mesh.TETS(i,cols);
           MM(Mindex) = Ge(rows,cols);
       end
    end
    % 右侧向量
    for p = 1:nQuad
        % 计算积分坐标
        pxyz = lambda(p,1)*mesh.POS(mesh.TETS(i,1),:)...
            +  lambda(p,2)*mesh.POS(mesh.TETS(i,2),:)...
            +  lambda(p,3)*mesh.POS(mesh.TETS(i,3),:)...
            +  lambda(p,4)*mesh.POS(mesh.TETS(i,4),:);
        % 计算点pxyz处的电流密度
        Jp = [0,0,0];        
        if mesh.ELE_TAGS(base+i,2) == CoilTag
            % 计算角度
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

% 构造稀疏矩阵
% 对角线上的元素与上三角元素
diagAIndex = (rowAIndex == colAIndex);
upperAIndex = ~diagAIndex;
% 稀疏矩阵A的大小为Ne*Ne
A = sparse(rowAIndex(diagAIndex),colAIndex(diagAIndex),AA(diagAIndex),Ne,Ne);
AUpper = sparse(rowAIndex(upperAIndex),colAIndex(upperAIndex),AA(upperAIndex),Ne,Ne);
A = A + AUpper + AUpper';% 根据对称性补全大矩阵
M = sparse(rowMIndex,colMIndex,MM,Ne,Nn);
% 计算f
bt = bt.*repmat(volume,1,6);
f = accumarray(elem2edge(:),bt(:),[Ne 1]);

%% 施加边界条件
% 一般求解的时候，不需要特别的设置边界条件，因为有一个默认的几何边界，那么为了
% 设置边界，就需要找出这个边界。以前的思路是拿几何尺寸来判断，缺点就是你必须知道
% 边界的形状和尺寸，非常不方便。但是，再想一下，如果一个单元是边界，二维的时候，
% 这个棱边界只存在于一个分网单元中，那么就可以查找出现一次的边，那就是边界了。
% 三维的思路也是一样的，只不过要查找的是面单元。
% 生成四面体的所有面
allFace = [mesh.TETS(:,[2 4 3]);mesh.TETS(:,[1 3 4]);mesh.TETS(:,[1 4 2]);mesh.TETS(:,[1 2 3])];

Nfall = length(allFace);% 面的数目
% 不重复的face
[face, i2, j] = unique(sort(allFace,2),'rows','legacy');

i1(j(Nfall:-1:1)) = Nfall:-1:1;% 只有出现一次的编号才不会变 
i1 = i1';
bdFlag = zeros(Nfall,1,'uint8');
bdFaceidx = i1(i1==i2);% 边界上的face的索引
bdFlag(bdFaceidx) = 1;% Dirichlet边界条件
bdFlag = reshape(bdFlag,NT,4);% 转换为对应的单元形式
% 查找边界棱和节点
isBdEdge = false(Ne,1);
isBdNode = false(Nn,1);
isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
bdEdge = edge(isBdEdge,:);% 将二维数组转换为向量
isBdNode(bdEdge) = true;% 标记边界节点
% 未知的变量
freeEdge = find(~isBdEdge);
freeNode = find(~isBdNode);
% 设置Neumann边界条件

% 设置Dirichlet边界条件
u = zeros(Ne,1);
if ~isempty(bdEdge)
    u(isBdEdge) = 0;% 固定边界条件
%     f = f - (A - M)*u;
    f(isBdEdge) = u(isBdEdge);
end
g0 = -M'*u;
%% 求解
% reduce to free dofs
A  = A(freeEdge,freeEdge);
M  = M(freeEdge,freeNode);
f = f(freeEdge);
g0 = g0(freeNode);
Ni = size(freeNode,1);
bigA = [A M;...
        M' sparse(Ni,Ni)];
tmp = bigA\[f;g0];

%% 输出求解结果