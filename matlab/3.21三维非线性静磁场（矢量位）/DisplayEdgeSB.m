function [] = DisplayEdgeSB(i,BrickElementConnectivityTable,GlobalNodeNumbering)
% 20190828
% by Poofee
% 绘制四面体内部的棱单元形函数的分布

ElementNodeXYZ = GetElementGlobalCoordinatesXYZ(...
    i,...
    BrickElementConnectivityTable,...
    GlobalNodeNumbering);

X = ElementNodeXYZ(:,1);    % Local node global x coordinates
Y = ElementNodeXYZ(:,2);    % Local node global y coordinates
Z = ElementNodeXYZ(:,3);    % Local node global z coordinates

Lines = [1 2 3 4 4 1;...
    2 3 4 1 2 3];
line(X(Lines),Y(Lines),Z(Lines),'Color',[0 0 0]);
hold on;
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
W = zeros(6,3);
for gridi = 1:length(gridx)
    N(1) = TetraNodalBasis(1,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    N(2) = TetraNodalBasis(2,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    N(3) = TetraNodalBasis(3,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    N(4) = TetraNodalBasis(4,X,Y,Z,gridx(gridi),gridy(gridi),gridz(gridi));
    
    dN(1,:) = dTetraNodalBasis(1,X,Y,Z);
    dN(2,:) = dTetraNodalBasis(2,X,Y,Z);
    dN(3,:) = dTetraNodalBasis(3,X,Y,Z);
    dN(4,:) = dTetraNodalBasis(4,X,Y,Z);
    if abs(sum(abs(N))-1) < 1e-10
        A = N(1)*dN(2,:)-N(2)*dN(1,:);
        quiver3(gridx(gridi),gridy(gridi),gridz(gridi),A(1),A(2),A(3),0.000008);
    end
end

axis equal;

end