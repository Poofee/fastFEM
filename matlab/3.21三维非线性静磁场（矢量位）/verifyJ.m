% 20190830
% by Poofee
% 验证一下能不能对电流密度在单元内进行离散插值。
% 四面体 (0,0,0) (1,0,0) (0,1,0) (0,0,1)
close all;
element = [0,0,0;
           1,0,0;
           0,1,0;
           0,0,1];
% 绘制单元
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

s = 20;
[gridx,gridy,gridz] = meshgrid(linspace(xmin,xmax,s),...
    linspace(ymin,ymax,s),...
    linspace(zmin,zmax,s));
gridx = reshape(gridx,[numel(gridx),1]);
gridy = reshape(gridy,[numel(gridy),1]);
gridz = reshape(gridz,[numel(gridz),1]);

N = zeros(4,1);
dN = zeros(4,3);
% W = zeros(6,3);
J = [1;1;1;1;1;1];
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

