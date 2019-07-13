function [A,FixNL] = magsolve(mesh,deltat,Ak,step,FixNLk)
% ������
MAG_CIR = [1,3,4];% ��·
COIL = [2];% ��Ȧ
AIR = [7];% ����
MOBILE_AIR = [6];% �������˶��Ŀ���
CORE = [5];% ����
INFINITE = [8];% ����Զ����

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
for idomain=length(Domain):-1:1
    if Domain(idomain) == 6 || Domain(idomain) == 7 || Domain(idomain) == 8
        FixNL(idomain,:) = [];
    end
end
disp(['һ���� ',num2str(length(FixNL)),' ���̶�����'])

num_nodes = mesh.nbNod;
num_elements = mesh.nbTriangles;

CE = zeros(3,3);      % CE --- �����洢ÿ����Ԫ��ϵ������
D = zeros(3,3);
M = zeros(3,3);
S = zeros(num_nodes,num_nodes);%ȫ�־���
F1 = zeros(num_nodes,1);

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

CoilTurn = 225;
Scoil = 0.36e-3;
tau = CoilTurn/Scoil;
J(coildomain) = tau*24/3;%������Ȧ����ĵ����ܶȣ���������Ϊ0
Ys(coildomain) = 1/3;

mu0 = 4*pi*1e-7;%�����ŵ���
mu = mu0*ones(num_elements,1);%����ÿһ����Ԫ�Ĵŵ��ʣ���ʼ��Ϊ�����ŵ���

% ���ҷǱ߽��ϵĽڵ�
freenodes = find(abs(X)>1e-8 & sqrt(X.^2+Y.^2)<0.11-1e-3);
% plot(X,Y,'r*');
% hold on
% plot(X(freenodes),Y(freenodes),'bo');
% axis equal
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
tol = 1e-6;%�������

if step == 2
    return;
end

for count = 1:steps
    disp(['��ʼ�� ',num2str(count),' �������Ե���']);
    count_fix = 0;
    %     װ��
    for i=1:num_elements
        if find([1,3,4,5] == Domain(i))% ��������
            dvdB = getdvdB(B(i));
            sigma = 2e-7;
        else%��������
            dvdB = 0;
            sigma = 0;% ���������������಻����
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
                
                M(row,col) = sigma*AREA(i)*(1/12+(row==col)*1/12)/deltat;
                
                S(NL(i,row),NL(i,col)) = S(NL(i,row),NL(i,col)) + CE(row,col) + M(row,col);
                F1(NL(i,row)) = F1(NL(i,row)) + D(row,col)*A(NL(i,col));
            end
            
            if find(FIXED_MESH == Domain(i))
                % ��������һ����A��MӦ��Ϊ��������ģ���������Ϊ0�����Ժ�д
                F1(NL(i,row)) = F1(NL(i,row)) + sum(M(row,:).*Ak(FixNLk(count_fix,:))')/deltat;
                % ��Ȧ��Ӧ�������ⲿ������Ȧ����ģ��������������ĵ�����Ϊ��0����д��
                F1(NL(i,row)) = F1(NL(i,row)) - (tau*AREA(i)/3)^2*Ys(i)/deltat*Ak(FixNLk(count_fix,row));
                % ��������
                F1(NL(i,row)) = F1(NL(i,row)) + J(i)*AREA(i)/3;
            end
            
            S(NL(i,row),NL(i,row)) = S(NL(i,row),NL(i,row))-(tau*AREA(i)/3)^2*Ys(i)/deltat;
        end
    end
    A_old = A;
    A(freenodes) = S(freenodes,freenodes)\F1(freenodes);
    % �ж����
    error = norm((A_old - A))/norm(A);
    disp(['������ ',num2str(error)]);
    if error < tol
        break;
    end
    %     ��һ��������ʼ��
    S = S - S;
    F1 = F1 - F1;
    %����B
    Bx = sum(R.*A(NL),2)./AREA./ydot/2;
    By = sum(Q.*A(NL),2)./AREA./ydot/2;
    B = sqrt(Bx.*Bx + By.*By);
    %����mu��ֻ���·���������
    mu(DomainNL) = B(DomainNL)./arrayfun(@getH,B(DomainNL));
    
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
    title('A','FontName','Times New Roman','FontSize',15);
    set(gca,'FontName','Times New Roman','FontSize',15);
    set(gca,'FontName','Times New Roman','FontSize',15);
    contourf(qx,qy,qz,20);colorbar
    axis equal
    h = gcf;
    set(h,'position',[474 59 658 798]);
    drawnow
end

end