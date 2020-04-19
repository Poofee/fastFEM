function [] = tet_post_MAG_Magnetostatic_A(mesh,u,elem2edge,xmin,xmax,ymin,ymax,zmin,zmax,gridsize)
% 显示线圈区域的A的矢量分布，理论上是没有z分量
figure;
Aplot = gca;axis(Aplot,'equal');hold on;
title(Aplot,'线圈内磁势A的分布');
% 显示线圈区域的B的矢量分布
figure;
Bplot = gca;axis(Bplot,'equal');hold on;
title(Bplot,'线圈内磁感应强度B的分布');

% 绘制线圈附近的区域
xminBound = floor(xmin/gridsize);
xmaxBound = ceil(xmax/gridsize);
yminBound = floor(ymin/gridsize);
ymaxBound = ceil(ymax/gridsize);
zminBound = floor(zmin/gridsize);
zmaxBound = ceil(zmax/gridsize);
% 生成四面体内的网格
s = 3;
[gridxBound,gridyBound,gridzBound] = meshgrid(xminBound:xmaxBound,...
    yminBound:ymaxBound,...
    zminBound:zmaxBound);
sizeBound = size(gridxBound);
gridxBound = gridsize*reshape(gridxBound,[numel(gridxBound),1]);
gridyBound = gridsize*reshape(gridyBound,[numel(gridyBound),1]);
gridzBound = gridsize*reshape(gridzBound,[numel(gridzBound),1]);
Boutput = zeros(size(gridxBound,1),3);
Aoutput = zeros(size(gridxBound,1),3);
for i=1:mesh.nbTets
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(i,1:4),3);
    
    % 四面体的边界，整数化后
    tetxmin = floor(min(X)/gridsize);
    tetxmax = ceil(max(X)/gridsize);
    tetymin = floor(min(Y)/gridsize);
    tetymax = ceil(max(Y)/gridsize);
    tetzmin = floor(min(Z)/gridsize);
    tetzmax = ceil(max(Z)/gridsize);
    % 生成四面体内的网格
    s = 3;
    [gridx,gridy,gridz] = meshgrid(tetxmin:tetxmax,...
        tetymin:tetymax,...
        tetzmin:tetzmax);
    gridx = gridsize*reshape(gridx,[numel(gridx),1]);
    gridy = gridsize*reshape(gridy,[numel(gridy),1]);
    gridz = gridsize*reshape(gridz,[numel(gridz),1]);
    
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
    edgeVector = [X(Lines')*[-1;1],Y(Lines')*[-1;1],Z(Lines')*[-1;1]];
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
%     [XX, YY, ZZ ] = tetNedelec_XYZ( eleLen, si, xi, yi, zi);
%     rotA = tetNedelec_rot( D, XX, YY, ZZ, Ai);
%     dd = D / mu0 ;
    
    dN(1,:) = dTetraNodalBasis(1,X,Y,Z);
    dN(2,:) = dTetraNodalBasis(2,X,Y,Z);
    dN(3,:) = dTetraNodalBasis(3,X,Y,Z);
    dN(4,:) = dTetraNodalBasis(4,X,Y,Z);
    
    BB = 2*ones(1,6)*(Ai.*eleLen.*si*ones(1,3).*cross(dN(Lines(1,:),:),dN(Lines(2,:),:)));
    
    if norm([tetg(1),tetg(2)]) < 0.024 && abs(tetg(3))<0.023%mesh.ELE_TAGS(base+i,2) ~= CoilTag
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
                
%                 quiver3(Aplot,gridx(gridi),gridy(gridi),gridz(gridi),Aee(1),Aee(2),Aee(3),1/3);
%                 hold on;drawnow;
%                 axis(Aplot,'equal');
%                 
%                 quiver3(Bplot,gridx(gridi),gridy(gridi),gridz(gridi),BB(1),BB(2),BB(3),1/150);
%                 hold on;drawnow;
%                 axis(Bplot,'equal');
                
                % save data
                zIndex = round(gridz(gridi)/gridsize)-zminBound;
                yIndex = round(gridy(gridi)/gridsize)-yminBound;
                xIndex = round(gridx(gridi)/gridsize)-xminBound;
                outputIndex = zIndex*(sizeBound(1)*sizeBound(2))+xIndex*(sizeBound(2))+yIndex;
                Boutput(outputIndex,:) = BB;
                Aoutput(outputIndex,:) = Aee;
            end
        end
    end
end
% vtk output
makevtk_struc_grid(gridxBound,gridyBound,gridzBound,Aoutput(:,1),Aoutput(:,2),Aoutput(:,3),'Aoutput.vtk');
makevtk_struc_grid(gridxBound,gridyBound,gridzBound,Boutput(:,1),Boutput(:,2),Boutput(:,3),'Boutput.vtk');
end