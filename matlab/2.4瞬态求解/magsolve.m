function [A,FixNL,AREA,ik1,FixNLIndex] = magsolve(t,mesh,time,Ak,FixNLk,ik,AREA_0,FixNLIndexk,axResult,plotErrorA)
% ������
MAG_CIR = [1,3,4];% ��·
COIL = [2];% ��Ȧ
AIR = [7];% ����
MOBILE_AIR = [6];% �������˶��Ŀ���
CORE = [5];% ����
INFINITE = [8];% ����Զ����


timebeta = 0.5;% ʱ����ɢ��betaֵ
% ��������
COMPRESSIBLE_PART = [AIR,INFINITE];% ��ѹ������
FIXED_PART = [MAG_CIR,COIL,MOBILE_AIR];% �����ƶ�����
MOBILE_PART = [CORE];% ���ƶ�����
FIXED_MESH = [MAG_CIR,COIL,CORE];
NONLINEAR_PART = [MAG_CIR,CORE];

X = mesh.POS(:,1);
Y = mesh.POS(:,2);
NL = mesh.TRIANGLES(:,1:3);
Domain = mesh.ELE_TAGS((mesh.nbElm-mesh.nbTriangles+1):end,2);
FixNL = NL;
% ɾ����������
FixNLIndex = 1:1:length(Domain);
for idomain=length(Domain):-1:1
    if Domain(idomain) == 6 || Domain(idomain) == 7 || Domain(idomain) == 8
        FixNL(idomain,:) = [];
        FixNLIndex(idomain) = [];
    end
end
disp(['һ���� ',num2str(length(FixNL)),' ���̶�����'])

num_nodes = mesh.nbNod;
num_elements = mesh.nbTriangles;

CE = zeros(3,3);      % CE --- �����洢ÿ����Ԫ��ϵ������
D = zeros(3,3);
M = zeros(3,3);
bbT = zeros(3,3);
S = zeros(num_nodes,num_nodes);%ȫ�־���
Jacobi = zeros(num_nodes,num_nodes);
mT = zeros(num_nodes,num_nodes);
mS = zeros(num_nodes,num_nodes);
mD = zeros(num_nodes,1);
F1 = zeros(num_nodes,1);
B1 = zeros(num_nodes,1);

%AREA = zeros(num_elements,1);%���������
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
% ��֪��Ϊʲô���Ĭ�ϲ�������
AREA = 0.5 * abs(Q(:,1).*R(:,2) - Q(:,2).*R(:,1));%���������

%��ģ����ֻ��5������1��2Ϊ��������5Ϊ��Ȧ����3��4Ϊ��о��
%����һ���������⡣
J = zeros(num_elements,1);%��������ܶȾ���
Ys = zeros(num_elements,1);
coildomain = find(Domain == 2);%Ѱ����Ȧ����ĵ�Ԫ
% ע�� | ���ԶԾ�����в��� ||�򲻿��ԣ��ᱨ��
MAG_CIRdomain = find(Domain == 1 | Domain == 3 | Domain == 4);
COREdomain = find(Domain == 5);

CoilTurn = 225;
Scoil = 8e-3*45e-3;
tau = CoilTurn/Scoil;
J(coildomain) = ik(t-1);%������Ȧ����ĵ����ܶȣ���������Ϊ0
Ys(coildomain) = 1/3;

mu0 = 4*pi*1e-7;%�����ŵ���
mu = mu0*ones(num_elements,1);%����ÿһ����Ԫ�Ĵŵ��ʣ���ʼ��Ϊ�����ŵ���

% ���ҷǱ߽��ϵĽڵ�
% Ҳ����
totalEdge = int32([NL(:,[1 2]); NL(:,[1 3]); NL(:,[2 3]);]);
NtotalEdge = length(totalEdge);
% �ԱߵĽڵ��Ž��д����ӵ͵���
sortedTotalEdge = sort(totalEdge,2);
% edge��ȫ�ֵ��⣬��СNTx2�������ʼ��ͽ�����
% j������sortedTotalEdge�е�������edge�еı��
[edge,i2,j] = unique(sortedTotalEdge,'rows','legacy');
i1(j(NtotalEdge:-1:1)) = NtotalEdge:-1:1;
i1 = i1';
bdEdgeIndex = i1(i1 == i2);

isbdNode = true(num_nodes,1);
% ���ڱ�����⣬POS��������ĵ㲢��ȫ�Ƿ����ڵ�
% ע�ⲻҪ��gmsh�����нڵ㵱�����������нڵ㣬��Ϊ������ܰ�����һЩ���νڵ㡣
allNode = mesh.TRIANGLES(:,1:3);
allNode = allNode(:);
uniNode = unique(allNode);
isbdNode(uniNode) = false;% �ҵ�ȫ�������ڵ�
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
%% �����������ı߽�
bdEdgeCoil = findDomainBoundary(NL,coildomain);
bdEdgeMAG_CIR = findDomainBoundary(NL,MAG_CIRdomain);
bdEdgeCORE = findDomainBoundary(NL,COREdomain);

% ����̲���P59ҳ��y
ydot = zeros(num_elements,1);
for i=1:num_elements
    %����������������,ע��P59ҳ��ʽ�Ǵ���ģ�Ӧ��Ϊx
    if XL(i,1)+XL(i,2)<1e-10 || XL(i,2)+XL(i,3)<1e-10 || XL(i,1)+XL(i,3)<1e-10
        ydot(i) = mean(XL(i,:));
    else
        ydot(i) = 1.5/(1/(XL(i,1)+XL(i,2))+1/(XL(i,1)+XL(i,3))+1/(XL(i,2)+XL(i,3)));
    end
end
% �����Ե���
A = zeros(num_nodes,1);%ÿ���ڵ�Ĵ��ƣ���Գ�������
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
%% ���� fix  A���ٶ���һ������һ���Ĺ̶��������ǰ�����ͬ˳��洢��Ԫ�ġ�
% Ϊ�˱�֤�������ȷ�����һ�������򣬿����ǲ���˳��һ�¡�
% ��Ҫ˼·�ǱȽ������
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
tol = 1e-6;%�������

ik1 = ik(t-1);

% �ų����·����ϵ���
deltT = (time(t)-time(t-1));
for couple_iteration=1:steps
    % �ų����ֵĵ���
    disp(['��ʼ�� ',num2str(couple_iteration),' �������Ե���']);

    Jacobi = Jacobi - Jacobi;
    mT = mT - mT;
    mD = mD -mD;
    mS = mS - mS;
    %     װ��
    for i=1:num_elements
        if find([1,3,4,5] == Domain(i))% ��������
            dvdB = getdvdB(B(i));
            sigma = 1/(2e-7);
        else%��������
            dvdB = 0;
            sigma = 0;% ���������������಻����
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
                % ��������
                mD(NL(i,row)) = mD(NL(i,row)) + tau*AREA(i)/3;
            end
        end
    end
    % ��������
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
    %     ��һ��������ʼ��
    
    %����B
    Bx = sum(R.*A(NL),2)./AREA./ydot/2;
    By = sum(Q.*A(NL),2)./AREA./ydot/2;
    B = sqrt(Bx.*Bx + By.*By);
    %����mu��ֻ���·��������������ֵBHΪ0�Ļ������ܱ�0��
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
    % ���Ƹ�������ı߽�
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
    
    % �ж����
    errorA = norm((A_old - A))/norm(A);
    errorI = norm((ik1_old - ik1))/norm(ik1);
    
    if errorA < tol && errorI < tol
        break;
    end
    disp(['�� ',num2str(couple_iteration),' ������ֵ ',num2str(ik1),' A',' ��ͨ ']);%num2str(coilphi)]
    disp(['������ ',num2str(errorA),' ',num2str(errorI)]);
 
    xd = get(plotErrorA,'XData');
    xd = [xd,length(xd)+1];
    yd = get(plotErrorA,'YData');
    yd = [yd,errorA];
    set(plotErrorA,'XData',xd,'YData',yd);

    drawnow
end

end