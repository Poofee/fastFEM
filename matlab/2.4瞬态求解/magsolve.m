function [A,FixNL,AREA,ik1,FixNLIndex] = magsolve(t,mesh,time,Ak,FixNLk,ik,AREA_0,FixNLIndexk,axResult,plotErrorA)
% 区域编号
MAG_CIR = [1,3,4];% 磁路
COIL = [2];% 线圈
AIR = [7];% 空气
MOBILE_AIR = [6];% 随衔铁运动的空气
CORE = [5];% 衔铁
INFINITE = [8];% 无穷远区域


timebeta = 0.5;% 时间离散的beta值
% 物理区域
COMPRESSIBLE_PART = [AIR,INFINITE];% 可压缩区域
FIXED_PART = [MAG_CIR,COIL,MOBILE_AIR];% 不可移动区域
MOBILE_PART = [CORE];% 可移动区域
FIXED_MESH = [MAG_CIR,COIL,CORE];
NONLINEAR_PART = [MAG_CIR,CORE];

X = mesh.POS(:,1);
Y = mesh.POS(:,2);
NL = mesh.TRIANGLES(:,1:3);
Domain = mesh.ELE_TAGS((mesh.nbElm-mesh.nbTriangles+1):end,2);
FixNL = NL;
% 删掉空气部分
FixNLIndex = 1:1:length(Domain);
for idomain=length(Domain):-1:1
    if Domain(idomain) == 6 || Domain(idomain) == 7 || Domain(idomain) == 8
        FixNL(idomain,:) = [];
        FixNLIndex(idomain) = [];
    end
end
disp(['一共有 ',num2str(length(FixNL)),' 个固定网格'])

num_nodes = mesh.nbNod;
num_elements = mesh.nbTriangles;

CE = zeros(3,3);      % CE --- 用来存储每个单元的系数矩阵
D = zeros(3,3);
M = zeros(3,3);
bbT = zeros(3,3);
S = zeros(num_nodes,num_nodes);%全局矩阵
Jacobi = zeros(num_nodes,num_nodes);
mT = zeros(num_nodes,num_nodes);
mS = zeros(num_nodes,num_nodes);
mD = zeros(num_nodes,1);
F1 = zeros(num_nodes,1);
B1 = zeros(num_nodes,1);

%AREA = zeros(num_elements,1);%三角形面积
P = zeros(num_elements,3);
Q = zeros(num_elements,3);
R = zeros(num_elements,3);

XL = X(NL);
YL = Y(NL);

Q(:,1) = YL(:,2) - YL(:,3);
Q(:,2) = YL(:,3) - YL(:,1);
Q(:,3) = YL(:,1) - YL(:,2);

R(:,1) = XL(:,3) - XL(:,2);
R(:,2) = XL(:,1) - XL(:,3);
R(:,3) = XL(:,2) - XL(:,1);

P(:,1) = XL(:,2).*YL(:,3) - XL(:,3).*YL(:,2);
P(:,2) = XL(:,3).*YL(:,1) - XL(:,1).*YL(:,3);
P(:,3) = XL(:,1).*YL(:,2) - XL(:,2).*YL(:,1);
% 不知道为什么面积默认不是正的
AREA = 0.5 * abs(Q(:,1).*R(:,2) - Q(:,2).*R(:,1));%三角形面积

%该模型里只有5个域，域1，2为空气，域5为线圈，域3，4为铁芯。
%这是一个线性问题。
J = zeros(num_elements,1);%计算电流密度矩阵
Ys = zeros(num_elements,1);
coildomain = find(Domain == 2);%寻找线圈区域的单元
% 注意 | 可以对矩阵进行操作 ||则不可以，会报错
MAG_CIRdomain = find(Domain == 1 | Domain == 3 | Domain == 4);
COREdomain = find(Domain == 5);

CoilTurn = 225;
Scoil = 8e-3*45e-3;
tau = CoilTurn/Scoil;
J(coildomain) = ik(t-1);%设置线圈区域的电流密度，其他其余为0
Ys(coildomain) = 1/3;

mu0 = 4*pi*1e-7;%空气磁导率
mu = mu0*ones(num_elements,1);%保存每一个单元的磁导率，初始化为空气磁导率

% 查找非边界上的节点
% 也就是
totalEdge = int32([NL(:,[1 2]); NL(:,[1 3]); NL(:,[2 3]);]);
NtotalEdge = length(totalEdge);
% 对边的节点编号进行处理，从低到高
sortedTotalEdge = sort(totalEdge,2);
% edge，全局的棱，大小NTx2，棱的起始点和结束点
% j，保存sortedTotalEdge中的数据在edge中的编号
[edge,i2,j] = unique(sortedTotalEdge,'rows','legacy');
i1(j(NtotalEdge:-1:1)) = NtotalEdge:-1:1;
i1 = i1';
bdEdgeIndex = i1(i1 == i2);

isbdNode = true(num_nodes,1);
% 由于编号问题，POS变量保存的点并不全是分网节点
% 注意不要把gmsh的所有节点当做分网的所有节点，因为里面可能包含了一些几何节点。
allNode = mesh.TRIANGLES(:,1:3);
allNode = allNode(:);
uniNode = unique(allNode);
isbdNode(uniNode) = false;% 找到全部分网节点
isbdNode(totalEdge(bdEdgeIndex,:)) = true;

freenodes = find(~isbdNode);

bdEdge = totalEdge(bdEdgeIndex,:);
% for i=1:size(bdEdge,1)
%     line(mesh.POS(bdEdge(i,:),1),mesh.POS(bdEdge(i,:),2),'Color',[0 0 0],'LineStyle','-','Marker','*');
% end
% axis equal
% hold on
% plot(X,Y,'r*');
% hold on
% plot(X(isbdNode),Y(isbdNode),'bo');
% axis equal
%% 计算各个区域的边界
bdEdgeCoil = findDomainBoundary(NL,coildomain);
bdEdgeMAG_CIR = findDomainBoundary(NL,MAG_CIRdomain);
bdEdgeCORE = findDomainBoundary(NL,COREdomain);

% 计算教材上P59页的y
ydot = zeros(num_elements,1);
for i=1:num_elements
    %两个点在坐标轴上,注意P59页公式是错误的，应当为x
    if XL(i,1)+XL(i,2)<1e-10 || XL(i,2)+XL(i,3)<1e-10 || XL(i,1)+XL(i,3)<1e-10
        ydot(i) = mean(XL(i,:));
    else
        ydot(i) = 1.5/(1/(XL(i,1)+XL(i,2))+1/(XL(i,1)+XL(i,3))+1/(XL(i,2)+XL(i,3)));
    end
end
% 非线性迭代
A = zeros(num_nodes,1);%每个节点的磁势，轴对称情形下
AlastFix = zeros(num_nodes,1);
B = zeros(num_elements,1);
DomainNL = 1:length(Domain);
for idomainnl=length(Domain):-1:1
    if find(NONLINEAR_PART == Domain(idomainnl))
        ;
    else
        DomainNL(idomainnl) = [];
    end
end
%% 计算 fix  A，假定上一步和这一步的固定分网都是按照相同顺序存储单元的。
% 为了保证程序的正确，添加一个检测程序，看看是不是顺序一致。
% 主要思路是比较面积。
count_fix = 0;
for i=1:num_elements
    if find(FIXED_MESH == Domain(i))
        count_fix = count_fix + 1;
        if t > 2
            if abs(AREA(i)/AREA_0(FixNLIndexk(count_fix))-1) > 1e-10
               error('Mesh order Not right!'); 
            end
        end
        for row = 1:3
            AlastFix(NL(i,row)) = Ak(FixNLk(count_fix,row)); 
        end
    end
end

steps = 100;
tol = 1e-6;%收敛误差

ik1 = ik(t-1);

% 磁场与电路的耦合迭代
deltT = (time(t)-time(t-1));
for couple_iteration=1:steps
    % 磁场部分的迭代
    disp(['开始第 ',num2str(couple_iteration),' 步非线性迭代']);

    Jacobi = Jacobi - Jacobi;
    mT = mT - mT;
    mD = mD -mD;
    mS = mS - mS;
    %     装配
    for i=1:num_elements
        if find([1,3,4,5] == Domain(i))% 铁磁区域
            dvdB = getdvdB(B(i));
            sigma = 1/(2e-7);
        else%线性区域
            dvdB = 0;
            sigma = 0;% 除了铁磁区域其余不考虑
        end
        for row = 1:3
            for col = 1:3
                D(row,col) = dvdB/ydot(i)/ydot(i)/ydot(i)/AREA(i);
                D(row,col) = D(row,col)* sum((R(i,row)*R(i,:)+Q(i,row)*Q(i,:)).*A(NL(i,:))')/4/AREA(i);
                D(row,col) = D(row,col)* sum((R(i,col)*R(i,:)+Q(i,col)*Q(i,:)).*A(NL(i,:))')/4/AREA(i);
                CE(row,col) = (R(i,row)*R(i,col)+Q(i,row)*Q(i,col))/4/AREA(i)/mu(i)/ydot(i);
                
                M(row,col) = sigma/ydot(i)*AREA(i)*(1/12+(row==col)*1/12);
                
                Jacobi(NL(i,row),NL(i,col)) = Jacobi(NL(i,row),NL(i,col)) + D(row,col) + CE(row,col);
                mT(NL(i,row),NL(i,col)) = mT(NL(i,row),NL(i,col)) + M(row,col);
                mS(NL(i,row),NL(i,col)) = mS(NL(i,row),NL(i,col)) + CE(row,col);
            end
            
            if Domain(i) == 2
                % 传导电流
                mD(NL(i,row)) = mD(NL(i,row)) + tau*AREA(i)/3;
            end
        end
    end
    % 计算大矩阵
    S = [timebeta*Jacobi+mT/deltT, -timebeta*mD;
        -timebeta*mD', -deltT*timebeta^2/2/pi*3];
    S1 = -[timebeta*mS+mT/deltT, -timebeta*mD;
        -timebeta*mD', -deltT*timebeta^2/2/pi*3;];
    S2 = [(timebeta-1)*mS+mT/deltT, -(timebeta-1)*mD;
        -timebeta*mD', -deltT*timebeta/2/pi*(timebeta-1)*3];
    F2 = S1*[A;ik1]+S2*[AlastFix;ik(t-1)]+[zeros(num_nodes,1);-deltT*timebeta/2/pi*24];
    A_old = A;
    ik1_old = ik1;
    XX = S([freenodes;end],[freenodes;end])\F2([freenodes;end]);  
    A(freenodes) = XX(1:end-1)+A_old(freenodes);
    ik1 = XX(end)+ik1_old;
    %     下一步迭代初始化
    
    %更新B
    Bx = sum(R.*A(NL),2)./AREA./ydot/2;
    By = sum(Q.*A(NL),2)./AREA./ydot/2;
    B = sqrt(Bx.*Bx + By.*By);
    %更新mu，只更新非线性区域，如果初值BH为0的话，可能被0除
    mu(DomainNL) = B(DomainNL)./arrayfun(@getH,B(DomainNL));
    
    FF = scatteredInterpolant(X,Y,A);
    tx = 0:1e-4:0.035;
    ty = -0.040:1e-4:0.050;
    [qx,qy] = meshgrid(tx,ty);
    qz = FF(qx,qy);
    
    cla(axResult);
    hold on
    title(axResult,['A,',num2str(time(t)),'s,',num2str(t),',',num2str(couple_iteration),','],'FontName','Times New Roman','FontSize',15);
    set(axResult,'FontName','Times New Roman','FontSize',15);
    set(axResult,'FontName','Times New Roman','FontSize',15);
    contourf(axResult,qx,qy,qz,20);colorbar(axResult);
%     axis equal
%     h = gcf;
%     sizeScreen = get(0,'ScreenSize');
%     width = sizeScreen(3);
%     height = sizeScreen(4);
%     set(h,'Position',[(width-0.8*height*0.7)/2 48 0.8*height*0.7 0.8*height]);
    % 绘制各个区域的边界
    for ie=1:size(bdEdgeCoil,1)
        line(axResult,mesh.POS(bdEdgeCoil(ie,:),1),mesh.POS(bdEdgeCoil(ie,:),2),'Color',[1 0 0],'LineStyle','-','Marker','.');
    end
    hold on;
    for ie=1:size(bdEdgeMAG_CIR,1)
        line(axResult,mesh.POS(bdEdgeMAG_CIR(ie,:),1),mesh.POS(bdEdgeMAG_CIR(ie,:),2),'Color',[1 0 0],'LineStyle','-','Marker','.');
    end
    for ie=1:size(bdEdgeCORE,1)
        line(axResult,mesh.POS(bdEdgeCORE(ie,:),1),mesh.POS(bdEdgeCORE(ie,:),2),'Color',[1 0 0],'LineStyle','-','Marker','.');
    end
%     axis equal
    drawnow
    
    % 判断误差
    errorA = norm((A_old - A))/norm(A);
    errorI = norm((ik1_old - ik1))/norm(ik1);
    
    if errorA < tol && errorI < tol
        break;
    end
    disp(['第 ',num2str(couple_iteration),' 步电流值 ',num2str(ik1),' A',' 磁通 ']);%num2str(coilphi)]
    disp(['迭代误差： ',num2str(errorA),' ',num2str(errorI)]);
 
    xd = get(plotErrorA,'XData');
    xd = [xd,length(xd)+1];
    yd = get(plotErrorA,'YData');
    yd = [yd,errorA];
    set(plotErrorA,'XData',xd,'YData',yd);

    drawnow
end

end