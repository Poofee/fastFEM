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
                   mesh.TETS(:,[2 3]); mesh.TETS(:,[4 2]); mesh.TETS(:,[3 4])]);
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
locEdge = [1 2;1 3;1 4;2 3;4 2;3 4;];
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
    % tetNedelec_mk
    [xi,yi,zi] = tetNedelec_mk(X,Y,Z);
    % tetNedelec_XYZ
    [XX, YY, ZZ ] = tetNedelec_XYZ( eleLen, si, xi, yi, zi);
    % tetNedelec_rotrot
    rr = tetNedelec_rotrot( D, XX, YY, ZZ) ;
    Ae = rr/mu0;
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
    Ge = tetNedelec_gradside( D, b, c, d, Fx, Fy, Fz);
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
        bt(i,:) = bt(i,:) + be_JoA';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%另一种计算方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 计算梯度
    gradN = zeros(4,3);
    gradN(1,:) = dTetraNodalBasis(1,X,Y,Z,0,0,0);
    gradN(2,:) = dTetraNodalBasis(2,X,Y,Z,0,0,0);
    gradN(3,:) = dTetraNodalBasis(3,X,Y,Z,0,0,0);
    gradN(4,:) = dTetraNodalBasis(4,X,Y,Z,0,0,0);
    % 计算梯度方法二：见 https://www.iue.tuwien.ac.at/phd/nentchev/node31.html
    gradN1 = zeros(4,3);
    gradN1(1,:) = cross(edgeVector(4,:),edgeVector(5,:))/D;
    gradN1(2,:) = cross(edgeVector(2,:),edgeVector(3,:))/D;
    gradN1(3,:) = cross(edgeVector(3,:),edgeVector(1,:))/D;
    gradN1(4,:) = cross(edgeVector(1,:),edgeVector(2,:))/D;
 
    % 计算 \nabla \times \mathbf{w}_{i}
    CurlW = RotWBasis(X,Y,Z,0,0,0).*(eleLen*ones(1,3));
    % 计算旋度方法二： https://www.iue.tuwien.ac.at/phd/nentchev/node41.html
    CurlW1 = eleLen*ones(1,3).*edgeVector(end:-1:1,:)*2/D1;
    % 棱单元矩阵，大小为6*6
    Ae1 = CurlW*CurlW'*volume(i)/mu0;
    % LM算子矩阵，大小为6*4
    Ge1 = (gradN(locEdge(:,2),:)-gradN(locEdge(:,1),:))*gradN(1:4,:)'*volume(i)/4.*(eleLen*ones(1,4));
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
            phi_k = phi_k * eleLen(k);
            rhs = dot(phi_k,Jp);
            bt1(i,k) = bt1(i,k) + w(p)*rhs*volume(i);
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
isBdNode = true(Nn,1);
% 由于编号问题，POS变量保存的点并不全是分网节点
allNode = mesh.TETS(:,1:4);
allNode = allNode(:);
uniNode = unique(allNode);
isBdNode(uniNode) = false;% 找到全部分网节点
% 将边界面上的点设为边界点
isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
bdEdge = edge(isBdEdge,:);
isBdNode(bdEdge) = true;% 标记边界节点
% 绘制边界，验证
hold on
for i=1:size(bdEdge,1)
    line(mesh.POS(bdEdge(i,:),1),mesh.POS(bdEdge(i,:),2),mesh.POS(bdEdge(i,:),3),'Color',[0 0 0]);
end
DisplayNodes(mesh.POS(isBdNode,:));
view(45,atan(1/sqrt(2))*180/pi);
% 未知的变量
freeEdge = find(~isBdEdge);
freeNode = find(~isBdNode);
fprintf('自由棱数目为%d,自由节点数目为%d\n',size(freeEdge,1),size(freeNode,1));
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
% 显示线圈区域的A的矢量分布，理论上是没有z分量
figure;
Aplot = gca;
title('线圈内磁势A的分布');
axis(Aplot,'equal');
% 显示线圈区域的B的矢量分布
figure;
Bplot = gca;axis(Bplot,'equal');
title('线圈内磁感应强度B的分布');
for i=1:mesh.nbTets
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(i,1:4),3);
    
    % 四面体的边界
    xmin = min(X);
    xmax = max(X);
    ymin = min(Y);
    ymax = max(Y);
    zmin = min(Z);
    zmax = max(Z);
    % 生成四面体内的网格
    s = 3;
    [gridx,gridy,gridz] = meshgrid(linspace(xmin,xmax,s),...
        linspace(ymin,ymax,s),...
        linspace(zmin,zmax,s));
    gridx = reshape(gridx,[numel(gridx),1]);
    gridy = reshape(gridy,[numel(gridy),1]);
    gridz = reshape(gridz,[numel(gridz),1]);
    
    % si
    si = tetNedelec_Direction(X,Y,Z);
    % tet_simple_a
    a = tet_simple_a( X, Y, Z);
    % tet_simple_b
    b = tet_simple_b( X, Y, Z);
    % tet_simple_c
    c = tet_simple_c( X, Y, Z);
    % tet_simple_d
    d = tet_simple_d( X, Y, Z);
    
    
    Lines = [1 1 1 2 2 3;...
             2 3 4 3 4 4];
%     line(X(Lines),Y(Lines),Z(Lines),'Color',[0 0 0]);
    edgeVector = [X(locEdge)*[-1;1],Y(locEdge)*[-1;1],Z(locEdge)*[-1;1]];
    % 计算棱长    
    eleLen = vecnorm(edgeVector,2,2);
    lengthtet = mean(eleLen);
    hold on;
    % 获取6条棱的结果
    Ai = u(elem2edge(i,:));
    
    % 计算B
    % tet_Volume6
    D = tet_Volume6(X,Y,Z);
    % tet_Center
    tetg = tet_Center( X,Y,Z);
    % tetNedelec_XYZ
    [XX, YY, ZZ ] = tetNedelec_XYZ( eleLen, si, xi, yi, zi);
    rotA = tetNedelec_rot( D, XX, YY, ZZ, Ai);
    dd = D / mu0 ;
    
    dN(1,:) = dTetraNodalBasis(1,X,Y,Z);
    dN(2,:) = dTetraNodalBasis(2,X,Y,Z);
    dN(3,:) = dTetraNodalBasis(3,X,Y,Z);
    dN(4,:) = dTetraNodalBasis(4,X,Y,Z);
    
    BB = 2*ones(1,6)*(Ai.*eleLen.*si*ones(1,3).*cross(dN(Lines(1,:),:),dN(Lines(2,:),:)));
%     for ib=1:4
%         k = mesh.TETS(i,ib) ;
%         for jb=1:3
%             Bres(k,jb) = Bres(k,jb) + rotA(jb) * dd ;
%         end
%     end
    if norm([tetg(1),tetg(2)]) < 0.016
        quiver3(Bplot,tetg(1),tetg(2),tetg(3),BB(1),BB(2),BB(3),lengthtet/norm(BB)/2);
        drawnow
        hold on;
    end
    
    if mesh.ELE_TAGS(base+i,2) ~= CoilTag
        continue;
    end
    maxerror = 0;
    % 绘制网格中每个点的箭头
    for gridi = 1:length(gridx)
        N(1) = TetraNodalBasis(1,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
        N(2) = TetraNodalBasis(2,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
        N(3) = TetraNodalBasis(3,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
        N(4) = TetraNodalBasis(4,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
        
        W = WBasis(X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi)).*(eleLen.*si*ones(1,3));
        
        % 判断点是否在四面体内
        if abs(sum(abs(N))-1) < 1e-10
            % 计算磁势A
            Aee = sum(W.*(Ai*ones(1,3)),1);
            quiver3(Aplot,gridx(gridi),gridy(gridi),gridz(gridi),Aee(1),Aee(2),Aee(3),lengthtet/norm(Aee)/s);
            drawnow
            hold on;
        end
    end
end
