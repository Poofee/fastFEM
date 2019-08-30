% 20190830
% by Poofee
% ��֤һ���ܲ��ܶԵ����ܶ��ڵ�Ԫ�ڽ�����ɢ��ֵ��
% ������ (0,0,0) (1,0,0) (0,1,0) (0,0,1)
% ���Ǽ򵥵�����������Ԫ���������������������ܴﵽ�Լ���Ҫ�������������Բ��ԽԶ��
% ÿһ��������ģԽ�ӽ�1����˹���ֵĵ����ƺ�����û��Ӱ�졣�Ҿ������������ʵ�����ˡ�
% ���ȡ���ڵ��λ�ã���Ԫ�Ĵ�С����posx=100ʱ������Ѿ�������ǧ��֮һ�ˡ����ԣ�����
% ���ǽ�������Ԫ���ֵ�ԽСԽ�ã��������z���ԽԶ�����ԽС���о�һ�µ�Ԫ�ĳ��ȣ����
% ������z��ľ���Ĵ��¹�ϵ���ҳ�һ����С�������ֵ��
close all;
posx = 1;
lengthtet = 1e-2;
element = [0+posx,0,0;
           lengthtet+posx,0,0;
           0+posx,lengthtet,0;
           0+posx,0,lengthtet];
% ���Ƶ�Ԫ
Lines = [1 2 3 4 4 1;...
    2 3 4 1 2 3];
X = element(:,1);
Y = element(:,2);
Z = element(:,3);
line(X(Lines),Y(Lines),Z(Lines),'Color',[0 0 0]);
axis equal
hold on
text(X,Y,Z,{'1','2','3','4'});

xmin = min(X);
xmax = max(X);
ymin = min(Y);
ymax = max(Y);
zmin = min(Z);
zmax = max(Z);

s = 30;
[gridx,gridy,gridz] = meshgrid(linspace(xmin,xmax,s),...
    linspace(ymin,ymax,s),...
    linspace(zmin,zmax,s));
gridx = reshape(gridx,[numel(gridx),1]);
gridy = reshape(gridy,[numel(gridy),1]);
gridz = reshape(gridz,[numel(gridz),1]);

N = zeros(4,1);
dN = zeros(4,3);
% W = zeros(6,3);
ip3 = [-0.7745966692414830,0.5555555555555550;
0.0000000000000000,0.8888888888888880;
0.7745966692414830,0.5555555555555550];
ip4 = [-0.8611363115940520,0.3478548451374530;
-0.3399810435848560,0.6521451548625460;
0.3399810435848560,0.6521451548625460;
0.8611363115940520,0.3478548451374530];

gausspoints = ip4(:,1);
weight = ip4(:,2);
J0 = 1;
J = zeros(6,1);
% ����ÿһ�����ϵ��
edgeMap = [1,2;2,3;3,1;1,4;2,4;3,4];
edgeMap = fliplr(edgeMap);
for i = 1:6
    x1 = element(edgeMap(i,1),1);
    y1 = element(edgeMap(i,1),2);
    z1 = element(edgeMap(i,1),3);
    x2 = element(edgeMap(i,2),1);
    y2 = element(edgeMap(i,2),2);
    z2 = element(edgeMap(i,2),3);
    xmiddle = 0.5*(x1+x2);
    ymiddle = 0.5*(y1+y2);
    zmiddle = 0.5*(z1+z2);
    
    x = xmiddle + (x2-x1)*0.5.*gausspoints;
    y = ymiddle + (y2-y1)*0.5.*gausspoints;
    z = zmiddle + (z2-z1)*0.5.*gausspoints;
    
    for ip=1:length(gausspoints)
        costheta(ip) = x(ip)/norm([x(ip),y(ip)]);
        if y(ip) > 0
            sintheta(ip) = sqrt(1-costheta(ip)*costheta(ip));
        else
            sintheta(ip) = -sqrt(1-costheta(ip)*costheta(ip));
        end
        J(i) = J(i) + J0 * (costheta(ip)*(y2-y1)*0.5-sintheta(ip)*(x2-x1)*0.5)*weight(ip);
    end
    
%     \int \vec{J} \cdot dl = \int \vec{J} \cdot (dx,dy)
%     =\int (-J sintheta,J costheta) \cdot ((x2-x1)*0.5*dt,(y2-y1)*0.5*dt)
%     =\int (-J sintheta * (x2-x1)*0.5 + J costheta * (y2-y1)*0.5) dt
%     -1<t<1
end
for gridi = 1:length(gridx)
    N(1) = TetraNodalBasis(1,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    N(2) = TetraNodalBasis(2,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    N(3) = TetraNodalBasis(3,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    N(4) = TetraNodalBasis(4,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    
    dN(1,:) = dTetraNodalBasis(1,X,Y,Z);
    dN(2,:) = dTetraNodalBasis(2,X,Y,Z);
    dN(3,:) = dTetraNodalBasis(3,X,Y,Z);
    dN(4,:) = dTetraNodalBasis(4,X,Y,Z);
    W = WBasis(X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    
    if abs(sum(abs(N))-1) < 1e-10
%         A = N(2)*dN(3,:)-N(3)*dN(2,:);
        A = sum(W.*(J*ones(1,3)),1);
        norm(A)
        quiver3(gridx(gridi),gridy(gridi),gridz(gridi),A(1),A(2),A(3),0.05);
    end
end

axis equal