function [A,FixNL,AREA,ik1] = magsolve(t,mesh,time,Ak,step,FixNLk,ik,phik)
% 区域编号
MAG_CIR = [1,3,4];% 磁路
COIL = [2];% 线圈
AIR = [7];% 空气
MOBILE_AIR = [6];% 随衔铁运动的空气
CORE = [5];% 衔铁
INFINITE = [8];% 无穷远区域

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
for idomain=length(Domain):-1:1
    if Domain(idomain) == 6 || Domain(idomain) == 7 || Domain(idomain) == 8
        FixNL(idomain,:) = [];
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
F1 = zeros(num_nodes,1);

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

CoilTurn = 225;
Scoil = 8e-3*45e-3;
tau = CoilTurn/Scoil;
J(coildomain) = ik(t-1);%设置线圈区域的电流密度，其他其余为0
Ys(coildomain) = 1/3;

mu0 = 4*pi*1e-7;%空气磁导率
mu = mu0*ones(num_elements,1);%保存每一个单元的磁导率，初始化为空气磁导率

% 查找非边界上的节点
freenodes = find(abs(X)>1e-8 & sqrt(X.^2+Y.^2)<0.11-1e-3);
% plot(X,Y,'r*');
% hold on
% plot(X(freenodes),Y(freenodes),'bo');
% axis equal
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
B = zeros(num_elements,1);
DomainNL = [1:length(Domain)];
for idomainnl=length(Domain):-1:1
    if find(NONLINEAR_PART == Domain(idomainnl))
        ;
    else
        DomainNL(idomainnl) = [];
    end
end
steps = 20;
tol = 1e-6;%收敛误差

ik1 = ik(t-1);

if step == 2
    return;
end
% 磁场与电路的耦合迭代

for couple_iteration=1:10000
    % 求解电路，得到新的电流
    % 计算这一步的线圈磁通
    coilphi = 0;
    for j=1:mesh.nbTriangles
        if(find(COIL == mesh.ELE_TAGS(mesh.nbElm-mesh.nbTriangles+j,2)))
            n = mesh.TRIANGLES(j,1:3);
            coilphi = coilphi + sum(A(n))*AREA(j)/3;
        end
    end
    coilphi = abs(coilphi * 2 * pi * tau);
    if couple_iteration == 1
        coilphi = phik;
    end
    iktmp = ik1;
    ik1 = 1/3*(24)+ 1/3*((phik - coilphi)/(time(t)-time(t-1)));
    alpha1 = 1e-1;
    if t > 16
        alpha1 = 1e-4;
    end
    if t > 28
        alpha1 = 1e-5;
    end
    ik1 = iktmp - alpha1*(iktmp-ik1);
    disp(['第 ',num2str(couple_iteration),' 步电流值 ',num2str(ik1),' A',' 磁通 ',num2str(coilphi)]);
    J(coildomain) = ik1;
    ik(t) = ik1;
    if abs(ik1 - iktmp)/iktmp < 1e-6
        break;
    end
    % 磁场部分的迭代
    for count = 1:steps
        disp(['开始第 ',num2str(count),' 步非线性迭代']);
        count_fix = 0;
        %     装配
        for i=1:num_elements
            if find([1,3,4,5] == Domain(i))% 铁磁区域
                dvdB = getdvdB(B(i));
                sigma = 1/2e-7;
            else%线性区域
                dvdB = 0;
                sigma = 0;% 除了铁磁区域其余不考虑
            end
            if find(FIXED_MESH == Domain(i))
                count_fix = count_fix + 1;
            end
            for row = 1:3
                for col = 1:3
                    D(row,col) = dvdB/ydot(i)/ydot(i)/ydot(i)/AREA(i);
                    D(row,col) = D(row,col)* sum((R(i,row)*R(i,:)+Q(i,row)*Q(i,:)).*A(NL(i,:))')/4/AREA(i);
                    D(row,col) = D(row,col)* sum((R(i,col)*R(i,:)+Q(i,col)*Q(i,:)).*A(NL(i,:))')/4/AREA(i);
                    CE(row,col) = D(row,col) + (R(i,row)*R(i,col)+Q(i,row)*Q(i,col))/4/AREA(i)/mu(i)/ydot(i);
                    
                    M(row,col) = sigma/ydot(i)*AREA(i)*(1/12+(row==col)*1/12)/(time(t)-time(t-1));
%                     bbT(row,col) = (tau*AREA(i)/3)^2*Ys(i)/deltat;
                    
                    S(NL(i,row),NL(i,col)) = S(NL(i,row),NL(i,col)) + CE(row,col) + M(row,col);
                    F1(NL(i,row)) = F1(NL(i,row)) + D(row,col)*A(NL(i,col));
                end
                
                if find(FIXED_MESH == Domain(i))
                    % 涡流，上一步的A，M应当为铁磁区域的，其余区域为0，可以合写
                    F1(NL(i,row)) = F1(NL(i,row)) + sum(M(row,:).*Ak(FixNLk(count_fix,:))');
                    % 线圈感应电流，这部分是线圈区域的，这里把其他区域的电阻设为了0，合写了
%                     F1(NL(i,row)) = F1(NL(i,row)) - sum(bbT(row,:).*Ak(FixNLk(count_fix,:))');
                    % 传导电流
                    F1(NL(i,row)) = F1(NL(i,row)) + J(i)*tau*AREA(i)/3;
                end
            end
        end
        A_old = A;
        A(freenodes) = S(freenodes,freenodes)\F1(freenodes);
        
        %     下一步迭代初始化
        S = S - S;
        F1 = F1 - F1;
        %更新B
        Bx = sum(R.*A(NL),2)./AREA./ydot/2;
        By = sum(Q.*A(NL),2)./AREA./ydot/2;
        B = sqrt(Bx.*Bx + By.*By);
        %更新mu，只更新非线性区域
        mu(DomainNL) = B(DomainNL)./arrayfun(@getH,B(DomainNL));
        [MM,II] = min(mu(DomainNL));
%         disp(MM);
        FF = scatteredInterpolant(X,Y,A);
        % figure
        % F = scatteredInterpolant(X,Y,A(1:num_nodes));
        tx = 0:1e-4:0.035;
        ty = -0.040:1e-4:0.050;
        [qx,qy] = meshgrid(tx,ty);
        qz = FF(qx,qy);
        % mesh(qx,qy,qz);
        % surf(qx,qy,qz);
        clf
        hold on
        title(['A,',num2str(time(t)),'s,',num2str(t),',',num2str(couple_iteration),',',num2str(count)],'FontName','Times New Roman','FontSize',15);
        set(gca,'FontName','Times New Roman','FontSize',15);
        set(gca,'FontName','Times New Roman','FontSize',15);
        contourf(qx,qy,qz,20);colorbar
        axis equal
        h = gcf;
        size = get(0,'ScreenSize');
        width = size(3);
        height = size(4);
        set(h,'Position',[(width-0.8*height*0.7)/2 48 0.8*height*0.7 0.8*height]);
        drawnow
        
        % 判断误差
        error = norm((A_old - A))/norm(A);
%         disp(['迭代误差： ',num2str(error)]);
        if error < tol
            break;
        end
    end
    
end


end