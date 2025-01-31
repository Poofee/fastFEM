% 注意：无法找到函数时请运行setpath
% 20200320 by Poofee
% 采用矢量棱单元法对三维静磁场进行求解
% 求解模型为线圈产生的磁场
% 参考adventure的代码进行编写
% 线圈的一些参数：
% 线圈电压：25 V
% 线圈匝数：1000
% 线圈内径：16 mm
% 线圈外径：24 mm
% 线圈高度：46 mm
% 20191207 by Poofee
% 调试寻找原因，目前就是算出来的LM非常大，A非常小，不是很合理。
% 那是不是f算错了？

% 想要调用tet的函数需要将路径加入path，或运行setpath
addpath(genpath('./../libtet'),'-begin');

close all
clear all
t0 = cputime;
%% 调用gmsh分网
fprintf('开始分网......\n');
% cmd = 'gmsh.exe  -3 -format msh2 coil.geo';
% [status,cmdout] = system(cmd);
% if status == 1
%     fprintf('未找到gmsh，请将gmsh.exe添加到PATH或者放置到当前目录.\n');
% end
mesh = load_gmsh2('coil.msh');
mesht = cputime;
fprintf('分网结束. 共 %d 个单元. 包括 %d 个节点, %d 个棱, %d 个三角形, %d 个四面体，用时 %.2f 秒\n',...
    mesh.nbElm,mesh.nbPoints,mesh.nbLines,mesh.nbTriangles,mesh.nbTets,mesht-t0);

%% 读取COMSOL分网
% mesh = ReadCOMSOL3D('coil.mphtxt');
%% 绘制
% DisplayNodes(mesh.POS);

% DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);

% DisplayEdgeSB(100,mesh.TETS(:,1:4),mesh.POS);

%% 物理参数
disp('设置参数...');
Ucoil = 20;%线圈电压
Ncoil = 3000;%线圈匝数
Rc = 1;%线圈电阻
Scoil = 0.046 * (0.024 - 0.016);%线圈面积
mu = 1;
mu0 = 4*pi*1e-7;%空气磁导率
Js = Ncoil*Ucoil/Rc/Scoil;%线圈电流密度，乘上mu0是为了方便计算
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
elem2edgeSign1 = zeros(NT,6);

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

p_elem = 4;
mp_elem  = 6;
nd_elem  = 10;
dimension = 3;

Bres = zeros(Nn,3);

ae = zeros(6,6);
be = zeros(6,1);
Joe = zeros(4,3);
maxabsdiffbt = zeros(mesh.nbTets,1);
maxabsdiffAe = zeros(mesh.nbTets,1);
% 计算右侧向量
% 高斯积分
[lambda,w] = quadpts3(2);
nQuad = size(lambda,1);
% 编号的顺序好像很重要
locEdge = [1 2;1 3;1 4;2 3;2 4;3 4;];
bt = zeros(NT,6);
bt1 = zeros(NT,6);

tmpIndex = 0;
mesh.POS = mesh.POS/1000;% gmsh的数据单位似乎是cm
fprintf('开始单元矩阵计算...\n');
for i=1:mesh.nbTets
    % tet_pickup_coordinate_4vertex
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(i,1:4),3);
    
    % 计算单元体积 tet_Volume6
    tmp = [ones(4,1),X,Y,Z];
    volume(i) = det(tmp)/6;
    if volume(i) < 0
        volume(i) = - volume(i);
        fprintf('警告：单元%d编号不规范\n',i);
    end
    edgeVector = [X(locEdge)*[-1;1],Y(locEdge)*[-1;1],Z(locEdge)*[-1;1]];
    % 另外一种计算体积的方法，利用向量
    D1 = dot(cross(edgeVector(1,:),edgeVector(2,:)),edgeVector(3,:));
    % tet_Volume6
    D = tet_Volume6(X,Y,Z);
    % 计算棱长 tet_SideLength
    eleLen = vecnorm(edgeVector,2,2);
    % tetNedelec_Direction
    si = tetNedelec_Direction(X,Y,Z);
    elem2edgeSign1(i,:) = si';
%     si = ones(6,1);
    % tetNedelec_mk
    [xi,yi,zi] = tetNedelec_mk(X,Y,Z);
    % tetNedelec_XYZ
    [XX, YY, ZZ ] = tetNedelec_XYZ( eleLen, si, xi, yi, zi);
    % tetNedelec_rotrot
    rr = tetNedelec_rotrot( D, XX, YY, ZZ) ;
    Ae1 = rr/mu0;
    % tetNedelec_km_mk
    [xy, yz, zx ] = tetNedelec_km_mk( X,Y,Z);
    % tet_Center
    tetg = tet_Center( X,Y,Z);
    % tetNedelec_Fxyz
    [Fx, Fy, Fz] = tetNedelec_Fxyz( eleLen, si, tetg, xi, yi, zi, xy, yz, zx);
    % tet_simple_a
    a = tet_simple_a( X, Y, Z);
    % tet_simple_b
    b = tet_simple_b( X, Y, Z);
    % tet_simple_c
    c = tet_simple_c( X, Y, Z);
    % tet_simple_d
    d = tet_simple_d( X, Y, Z) ;
    % tetNedelec_gradside
    Ge1 = tetNedelec_gradside( D, b, c, d, Fx, Fy, Fz);
    % 使用adventure Newton-Cotes method的代码计算b
    if mesh.ELE_TAGS(base+i,2) == CoilTag
        
        % 计算四个顶点的电流密度
        for iele=1:4
            pxyz = mesh.POS(mesh.TETS(i,iele),:);
            % 计算角度
            costheta = pxyz(1)/norm(pxyz(1:2));
            sintheta = sqrt(1-costheta*costheta);
            if pxyz(2) < 0
                sintheta = -sintheta;
            end
            Joe(iele,:) = Js*[-sintheta,costheta,0];
        end
        be_JoA = tetNedelec_JoA( Joe, si, eleLen, b, c, d);
        bt1(i,:) = bt1(i,:) + be_JoA';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%另一种计算方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 计算梯度
    gradN = dTetraNodalBasis(X,Y,Z,0,0,0);
    % 计算梯度方法二：见 https://www.iue.tuwien.ac.at/phd/nentchev/node31.html
    gradN1 = zeros(4,3);
    gradN1(1,:) = cross(edgeVector(4,:),edgeVector(5,:))/D;
    gradN1(2,:) = cross(edgeVector(2,:),edgeVector(3,:))/D;
    gradN1(3,:) = cross(edgeVector(3,:),edgeVector(1,:))/D;
    gradN1(4,:) = cross(edgeVector(1,:),edgeVector(2,:))/D;
    
    % 计算 \nabla \times \mathbf{w}_{i}
    CurlW = RotWBasis(X,Y,Z,0,0,0).*(eleLen.*elem2edgeSign(i,:)'*ones(1,3));
    % 计算旋度方法二： https://www.iue.tuwien.ac.at/phd/nentchev/node41.html
    CurlW1 = eleLen.*elem2edgeSign(i,:)'*ones(1,3).*edgeVector(end:-1:1,:)*2/D1;
    % 棱单元矩阵，大小为6*6
    Ae = CurlW*CurlW'*volume(i)/mu0;
    % LM算子矩阵，大小为6*4
    Ge = (gradN(locEdge(:,2),:)-gradN(locEdge(:,1),:))*gradN(1:4,:)'*volume(i)/4.*(eleLen.*elem2edgeSign(i,:)'*ones(1,4));
    % 保存A的位置矩阵
    for rows=1:6
        for cols=rows:6
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
    
    % 使用gauess积分法求解右侧向量
    for p = 1:nQuad
        % 计算积分坐标
        pxyz = lambda(p,:)*mesh.POS(mesh.TETS(i,1:4),:);
        % 计算点pxyz处的电流密度
        Jp = [0,0,0];
        if mesh.ELE_TAGS(base+i,2) == CoilTag
            % 计算角度
            costheta = pxyz(1)/norm(pxyz(1:2));
            sintheta = sqrt(1-costheta*costheta);
            if pxyz(2) < 0
                sintheta = -sintheta;
            end
            Jp = Js*[-sintheta,costheta,0];
        end
        %         figure(2);
        %         quiver3(pxyz(1),pxyz(2),pxyz(3),Jp(1),Jp(2),Jp(3),1e-2);
        %         axis equal
        %         drawnow
        %         hold on;
        for k=1:6
            a1 = locEdge(k,1);
            a2 = locEdge(k,2);
            phi_k = lambda(p,a1)*gradN(a2,:)-lambda(p,a2)*gradN(a1,:);
            phi_k = phi_k * eleLen(k)*elem2edgeSign(i,k);
            rhs = dot(phi_k,Jp);
            bt(i,k) = bt(i,k) + w(p)*rhs*volume(i);
        end
    end
    diffbt = (abs(bt1(i,:)) - abs(bt(i,:)))./abs(bt(i,:));
    maxabsdiffbt(i) = max(abs(diffbt));
    diffAe = (abs(Ae) - abs(Ae1))./abs(Ae);
    maxabsdiffAe(i) = max(abs(diffAe),[],'all');
    
    if max(abs(diffbt))> 5
        disp('too big');
    end
    if maxabsdiffAe(i)> 5
        disp('too big');
    end
end
fprintf('开始进行矩阵装配...\n');
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
% bt = bt.*repmat(volume,1,6);
f = accumarray(elem2edge(:),bt(:),[Ne 1]);
f1 = accumarray(elem2edge(:),bt1(:),[Ne 1]);
%% 施加边界条件
fprintf('开始设置边界条件...\n');
[isBdEdge,isBdNode,freeEdge,freeNode] = findDomain3DBoundary(mesh.TETS,edge,elem2edge);
% 设置Neumann边界条件

% 设置Dirichlet边界条件
u = zeros(Ne,1);
if ~isempty(edge(isBdEdge,:))
    u(isBdEdge) = 0;% 固定边界条件
    %     f = f - (A - M)*u;
    f(isBdEdge) = u(isBdEdge);
end
g0 = -M'*u;
%% 求解
fprintf('开始求解...\n');
% reduce to free dofs
A  = A(freeEdge,freeEdge);
M  = M(freeEdge,freeNode);
f = f(freeEdge);
f1 = f1(freeEdge);
g0 = g0(freeNode);
Ni = size(freeNode,1);
Nei = size(freeEdge,1);
bigA = [A M;...
    M' sparse(Ni,Ni)];
fprintf('矩阵bigA的秩为%d,条件数为%d,大小为%d\n',sprank(bigA),condest(bigA),size(bigA,1));
sols = bigA\[f;g0];
% 检查误差
p = zeros(Nn,1);
% u(freeEdge) = bicg(A,f,1e-6,1000);% 磁场的解
u(freeEdge) = sols(1:Nei);% 磁场的解
p(freeNode) = sols(Nei+1:Nei+Ni);
residual = norm([f;g0] - bigA*sols);
fprintf('反算误差为%f\n',residual);
% 检查Lagrange multiplier是否正确，正确范围是多少？
% LM算子是额外引入的，应该有M*p(freeNode)=0，这样才不会影响
% 原来的式子的结果。
normp = norm(p);
if(normp>1.0/Nn)
    fprintf('Lagrange multiplier大小为%f\n,计算结果不对\n',normp);
end
%% 输出求解结果
fprintf('绘制求解结果...\n');
% 设置网格的大小，因为网格是均匀的，所以应该能够严格的整数化
gridsize = 5e-3;
xmin = -24e-3;xmax = 24e-3;
ymin = -24e-3;ymax = 24e-3;
zmin = -23e-3;zmax = 23e-3;
tet_post_MAG_Magnetostatic_A(mesh,u,elem2edge,elem2edgeSign,xmin,xmax,ymin,ymax,zmin,zmax,gridsize);

