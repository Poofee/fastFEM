% 20200320 by Poofee
% ����ʸ���ⵥԪ������ά���ų��������
% ���ģ��Ϊ��Ȧ�����Ĵų�
% �ο�adventure�Ĵ�����б�д
% ��Ȧ��һЩ������
% ��Ȧ��ѹ��25 V
% ��Ȧ������1000
% ��Ȧ�ھ���16 mm
% ��Ȧ�⾶��24 mm
% ��Ȧ�߶ȣ�46 mm
% 20191207 by Poofee
% ����Ѱ��ԭ��Ŀǰ�����������LM�ǳ���A�ǳ�С�����Ǻܺ���
% ���ǲ���f����ˣ�
close all
clear all
t0 = cputime;
%% ����gmsh����
fprintf('��ʼ����......\n');
% cmd = 'gmsh.exe  -3 -format msh2 coil.geo';
% [status,cmdout] = system(cmd);
% if status == 1
%     fprintf('δ�ҵ�gmsh���뽫gmsh.exe��ӵ�PATH���߷��õ���ǰĿ¼.\n');
% end
mesh = load_gmsh2('coil.msh');
mesht = cputime;
fprintf('��������. �� %d ����Ԫ. ���� %d ���ڵ�, %d ����, %d ��������, %d �������壬��ʱ %.2f ��\n',...
    mesh.nbElm,mesh.nbPoints,mesh.nbLines,mesh.nbTriangles,mesh.nbTets,mesht-t0);

%% ��ȡCOMSOL����
% mesh = ReadCOMSOL3D('coil.mphtxt');
%% ����
% DisplayNodes(mesh.POS);

% DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);

% DisplayEdgeSB(100,mesh.TETS(:,1:4),mesh.POS);

%% �������
disp('���ò���...');
Ucoil = 20;%��Ȧ��ѹ
Ncoil = 3000;%��Ȧ����
Rc = 1;%��Ȧ����
Scoil = 0.046 * (0.024 - 0.016);%��Ȧ���
mu = 1;
mu0 = 4*pi*1e-7;%�����ŵ���
Js = Ncoil*Ucoil/Rc/Scoil;%��Ȧ�����ܶȣ�����mu0��Ϊ�˷������
AirTag = 80;
CoilTag = 1;

%% ��������Ŀ
NT = size(mesh.TETS,1);% ��Ԫ��Ŀ
Nn = size(mesh.POS,1);% �ڵ���Ŀ
% ���ɵ�Ԫ�����е��⣬˳���Ű���...
totalEdge = int32([mesh.TETS(:,[1 2]); mesh.TETS(:,[1 3]); mesh.TETS(:,[1 4]); ...
                   mesh.TETS(:,[2 3]); mesh.TETS(:,[4 2]); mesh.TETS(:,[3 4])]);
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
% �����Ҳ�����
% ��˹����
[lambda,w] = quadpts3(2);
nQuad = size(lambda,1);
% ��ŵ�˳��������Ҫ
locEdge = [1 2;1 3;1 4;2 3;4 2;3 4;];
bt = zeros(NT,6);
bt1 = zeros(NT,6);

tmpIndex = 0;
mesh.POS = mesh.POS/1000;% gmsh�����ݵ�λ�ƺ���cm
fprintf('��ʼ��Ԫ�������...\n');
for i=1:mesh.nbTets
    % tet_pickup_coordinate_4vertex
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(i,1:4),3);
    
    % ���㵥Ԫ��� tet_Volume6
    tmp = [ones(4,1),X,Y,Z];
    volume(i) = det(tmp)/6;
    if volume(i) < 0
       volume(i) = - volume(i);
       fprintf('���棺��Ԫ%d��Ų��淶\n',i);
    end
    edgeVector = [X(locEdge)*[-1;1],Y(locEdge)*[-1;1],Z(locEdge)*[-1;1]];
    % ����һ�ּ�������ķ������������� 
    D1 = dot(cross(edgeVector(1,:),edgeVector(2,:)),edgeVector(3,:));
    % tet_Volume6
    D = tet_Volume6(X,Y,Z);
    % �����ⳤ tet_SideLength    
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
    % ʹ��adventure Newton-Cotes method�Ĵ������b
    if mesh.ELE_TAGS(base+i,2) == CoilTag
        
        % �����ĸ�����ĵ����ܶ�
        for iele=1:4
            pxyz = mesh.POS(mesh.TETS(i,iele),:);
            % ����Ƕ�
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
    %%��һ�ּ��㷽��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % �����ݶ�
    gradN = zeros(4,3);
    gradN(1,:) = dTetraNodalBasis(1,X,Y,Z,0,0,0);
    gradN(2,:) = dTetraNodalBasis(2,X,Y,Z,0,0,0);
    gradN(3,:) = dTetraNodalBasis(3,X,Y,Z,0,0,0);
    gradN(4,:) = dTetraNodalBasis(4,X,Y,Z,0,0,0);
    % �����ݶȷ��������� https://www.iue.tuwien.ac.at/phd/nentchev/node31.html
    gradN1 = zeros(4,3);
    gradN1(1,:) = cross(edgeVector(4,:),edgeVector(5,:))/D;
    gradN1(2,:) = cross(edgeVector(2,:),edgeVector(3,:))/D;
    gradN1(3,:) = cross(edgeVector(3,:),edgeVector(1,:))/D;
    gradN1(4,:) = cross(edgeVector(1,:),edgeVector(2,:))/D;
 
    % ���� \nabla \times \mathbf{w}_{i}
    CurlW = RotWBasis(X,Y,Z,0,0,0).*(eleLen*ones(1,3));
    % �������ȷ������� https://www.iue.tuwien.ac.at/phd/nentchev/node41.html
    CurlW1 = eleLen*ones(1,3).*edgeVector(end:-1:1,:)*2/D1;
    % �ⵥԪ���󣬴�СΪ6*6
    Ae1 = CurlW*CurlW'*volume(i)/mu0;
    % LM���Ӿ��󣬴�СΪ6*4
    Ge1 = (gradN(locEdge(:,2),:)-gradN(locEdge(:,1),:))*gradN(1:4,:)'*volume(i)/4.*(eleLen*ones(1,4));
    % ����A��λ�þ���
    for rows=1:6
       for cols=rows:6
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
    
    % ʹ��gauess���ַ�����Ҳ�����
    for p = 1:nQuad
        % �����������
        pxyz = lambda(p,:)*mesh.POS(mesh.TETS(i,1:4),:);
        % �����pxyz���ĵ����ܶ�
        Jp = [0,0,0];        
        if mesh.ELE_TAGS(base+i,2) == CoilTag
            % ����Ƕ�
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
fprintf('��ʼ���о���װ��...\n');
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
% bt = bt.*repmat(volume,1,6);
f = accumarray(elem2edge(:),bt(:),[Ne 1]);
f1 = accumarray(elem2edge(:),bt1(:),[Ne 1]);
%% ʩ�ӱ߽�����
fprintf('��ʼ���ñ߽�����...\n');
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
isBdNode = true(Nn,1);
% ���ڱ�����⣬POS��������ĵ㲢��ȫ�Ƿ����ڵ�
allNode = mesh.TETS(:,1:4);
allNode = allNode(:);
uniNode = unique(allNode);
isBdNode(uniNode) = false;% �ҵ�ȫ�������ڵ�
% ���߽����ϵĵ���Ϊ�߽��
isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
bdEdge = edge(isBdEdge,:);
isBdNode(bdEdge) = true;% ��Ǳ߽�ڵ�
% ���Ʊ߽磬��֤
hold on
for i=1:size(bdEdge,1)
    line(mesh.POS(bdEdge(i,:),1),mesh.POS(bdEdge(i,:),2),mesh.POS(bdEdge(i,:),3),'Color',[0 0 0]);
end
DisplayNodes(mesh.POS(isBdNode,:));
view(45,atan(1/sqrt(2))*180/pi);
% δ֪�ı���
freeEdge = find(~isBdEdge);
freeNode = find(~isBdNode);
fprintf('��������ĿΪ%d,���ɽڵ���ĿΪ%d\n',size(freeEdge,1),size(freeNode,1));
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
fprintf('��ʼ���...\n');
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
fprintf('����bigA����Ϊ%d,������Ϊ%d,��СΪ%d\n',sprank(bigA),condest(bigA),size(bigA,1));
sols = bigA\[f;g0];
% ������
p = zeros(Nn,1);
% u(freeEdge) = bicg(A,f,1e-6,1000);% �ų��Ľ�
u(freeEdge) = sols(1:Nei);% �ų��Ľ�
p(freeNode) = sols(Nei+1:Nei+Ni);
residual = norm([f;g0] - bigA*sols);
fprintf('�������Ϊ%f\n',residual);
% ���Lagrange multiplier�Ƿ���ȷ����ȷ��Χ�Ƕ��٣�
% LM�����Ƕ�������ģ�Ӧ����M*p(freeNode)=0�������Ų���Ӱ��
% ԭ����ʽ�ӵĽ����                
normp = norm(p);
if(normp>1.0/Nn)
    fprintf('Lagrange multiplier��СΪ%f\n,����������\n',normp);
end
%% ��������
fprintf('���������...\n');
% ��ʾ��Ȧ�����A��ʸ���ֲ�����������û��z����
figure;
Aplot = gca;
title('��Ȧ�ڴ���A�ķֲ�');
axis(Aplot,'equal');
% ��ʾ��Ȧ�����B��ʸ���ֲ�
figure;
Bplot = gca;axis(Bplot,'equal');
title('��Ȧ�ڴŸ�Ӧǿ��B�ķֲ�');
for i=1:mesh.nbTets
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(i,1:4),3);
    
    % ������ı߽�
    xmin = min(X);
    xmax = max(X);
    ymin = min(Y);
    ymax = max(Y);
    zmin = min(Z);
    zmax = max(Z);
    % �����������ڵ�����
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
    % �����ⳤ    
    eleLen = vecnorm(edgeVector,2,2);
    lengthtet = mean(eleLen);
    hold on;
    % ��ȡ6����Ľ��
    Ai = u(elem2edge(i,:));
    
    % ����B
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
    % ����������ÿ����ļ�ͷ
    for gridi = 1:length(gridx)
        N(1) = TetraNodalBasis(1,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
        N(2) = TetraNodalBasis(2,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
        N(3) = TetraNodalBasis(3,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
        N(4) = TetraNodalBasis(4,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
        
        W = WBasis(X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi)).*(eleLen.*si*ones(1,3));
        
        % �жϵ��Ƿ�����������
        if abs(sum(abs(N))-1) < 1e-10
            % �������A
            Aee = sum(W.*(Ai*ones(1,3)),1);
            quiver3(Aplot,gridx(gridi),gridy(gridi),gridz(gridi),Aee(1),Aee(2),Aee(3),lengthtet/norm(Aee)/s);
            drawnow
            hold on;
        end
    end
end
