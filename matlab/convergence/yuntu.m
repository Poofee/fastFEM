% clear all
% close all
%2017-10-07
%实现四边形求解后的云图的绘制
%plot the figure


PIE = 4.0*atan(1);      % PIE --- 圆周率π
U0 = 4*PIE*1e-7;   % U0 --- 真空磁导率

clear all
close all

dims = 2;
nnodes = 4;
% 获得单元和节点关系信息，并且得到每个几何Domain与其单元的关系
fp = fopen('D:\fastFEM\model\35928quad.mphtxt', 'r');
%-------------Read the head
for i=1:1:18
    tline = fgets(fp);
end
%--------------mesh point
num_nodes = fscanf(fp, '%d # number of mesh points\n', [1,1]);
pts_ind = fscanf(fp, '%d # lowest mesh point index\n', [1,1]);
fgets(fp);

xy = fscanf(fp, '%lf %lf \n', [dims,num_nodes]);
xy = xy';
%--------------vertex
for i=1:7
    fgets(fp);
end
num_vtx_ns = fscanf(fp, '%d # number of nodes per element\n', [1,1]);
num_vtx_ele = fscanf(fp, '%d # number of elements\n', [1,1]);
fgets(fp);
vtx = fscanf(fp, '%d \n', [1,num_vtx_ele]);
vtx = vtx';
%--------------vertex
num_vtx_ele2=fscanf(fp, '%d # number of geometric entity indices\n', [1,1]);
fgets(fp);
vtx2=fscanf(fp, '%d \n', [1,num_vtx_ele2]);
vtx2 = vtx2';
%--------------boundary
for i=1:5
    fgets(fp);
end
num_bdr_ns=fscanf(fp, '%d # number of nodes per element\n', [1,1]);
num_bdr_ele=fscanf(fp, '%d # number of elements\n', [1,1]);
fgets(fp);
p=fscanf(fp, '%d %d\n', [2,num_bdr_ele]);
p = p'+1;
NDPP=p(:,1);
%--------------entity
num_entity=fscanf(fp, '%d # number of geometric entity indices\n', [1,1]);
fgets(fp);
entity = fscanf(fp, '%d \n', [1,num_entity]);
entity = entity'+1;
%--------------elements
for i=1:5
    fgets(fp);
end
ns_per_ele=fscanf(fp, '%d # number of nodes per element\n', [1,1]);
num_elements=fscanf(fp, '%d # number of elements\n', [1,1]);
fgets(fp);
NL=fscanf(fp, '%d %d %d %d\n', [nnodes,num_elements]);
NL = NL';
NL=NL+1;
%--------------Domain
num_domain=fscanf(fp, '%d # number of geometric entity indices\n', [1,1]);
fgets(fp);
Domain = fscanf(fp, '%d \n', [1,num_domain]);
Domain = Domain';
fclose(fp);
 %COMSOL的读取格式貌似不是按照1234的顺序，而是
%按照1243的顺序，这里交换了一下。
t34 = NL(:,3);
NL(:,3) = NL(:,4);
NL(:,4) = t34;

Domain3 = find(Domain==3);
Domain4 = find(Domain>=4 & Domain<=7);
D34 = [Domain3;Domain4];

% 读入armadillo导出的数据
fp = fopen('NRpmA.txt','r');
fgets(fp);
Num = fscanf(fp,'%d %d',[1 2]);
A = fscanf(fp,'%lf',[Num(1) 1]);
fclose(fp);

FNR = scatteredInterpolant(xy(:,1),xy(:,2),A(1:num_nodes));
% figure
% F = scatteredInterpolant(X,Y,A(1:num_nodes));
tx = 0:0.001/4:0.025;
ty = -0.025:0.001/4:0.025;
[qx,qy] = meshgrid(tx,ty);
qz = FNR(qx,qy);
% mesh(qx,qy,qz);
% surf(qx,qy,qz);
figure
subplot(1,3,2);hold on
contourf(qx,qy,qz,20);
c = colorbar;
c.FontSize = 20;
title('NR','FontName','Times New Roman','FontSize',24);
ax = gca;
ax.FontName = 'Times New Roman';
xticklabels({});
yticklabels({});
for i=1:length(D34)
    fill(xy(NL(D34(i),:),1),xy(NL(D34(i),:),2),'w','EdgeAlpha',0,'FaceAlpha',0.3);
    hold on
%     plot([XL(D34(i),1:3),XL(D34(i),1)],[YL(D34(i),1:3),YL(D34(i),1)],'-g');
end
axis equal

% % 读入COMSOL的数据文件，进行云图绘制
% fp = fopen('comsolA.txt','r');
% % 读取前面没用的9行
% for i = 1:9
%     fgets(fp);
% end
% % 读取三列数据
% comsoldata = fscanf(fp,'%lf %lf %lf\n',[ 3 Num(1)]);
% comsoldata = comsoldata';
% fclose(fp);
% 
% FCOMSOL = scatteredInterpolant(comsoldata(:,1),comsoldata(:,2),comsoldata(:,3));
% % figure
% % F = scatteredInterpolant(X,Y,A(1:num_nodes));
% tx = 0:0.001/4:0.025;
% ty = -0.025:0.001/4:0.025;
% [qx,qy] = meshgrid(tx,ty);
% qz = FCOMSOL(qx,qy);
% % mesh(qx,qy,qz);
% % surf(qx,qy,qz);
% subplot(1,3,1);
% contourf(qx,qy,qz,20);
% c = colorbar;
% c.FontSize = 20;
% title('COMSOL','FontName','Times New Roman','FontSize',24);
% ax = gca;
% ax.FontName = 'Times New Roman';
% xticklabels({});
% yticklabels({});
% hold on
% for i=1:length(D34)
%     fill(xy(NL(D34(i),:),1),xy(NL(D34(i),:),2),'w','EdgeAlpha',0,'FaceAlpha',0.3);
%     hold on
% %     plot([XL(D34(i),1:3),XL(D34(i),1)],[YL(D34(i),1:3),YL(D34(i),1)],'-g');
% end
% axis equal

% 读入armadillo导出的数据
fp = fopen('TLMpmA.txt','r');
fgets(fp);
Num = fscanf(fp,'%d %d',[1 2]);
A = fscanf(fp,'%lf',[Num(1) 1]);
fclose(fp);

FTLM = scatteredInterpolant(xy(:,1),xy(:,2),A(1:num_nodes));
% figure
% F = scatteredInterpolant(X,Y,A(1:num_nodes));
tx = 0:0.001/4:0.025;
ty = -0.025:0.001/4:0.025;
[qx,qy] = meshgrid(tx,ty);
qz = FTLM(qx,qy);
% mesh(qx,qy,qz);
% surf(qx,qy,qz);
subplot(1,3,3);hold on
contourf(qx,qy,qz,20);
c = colorbar;
c.FontSize = 20;
title('TLM','FontName','Times New Roman','FontSize',24);
ax = gca;
ax.FontName = 'Times New Roman';
xticklabels({});
yticklabels({});
for i=1:length(D34)
    fill(xy(NL(D34(i),:),1),xy(NL(D34(i),:),2),'w','EdgeAlpha',0,'FaceAlpha',0.3);
    hold on
%     plot([XL(D34(i),1:3),XL(D34(i),1)],[YL(D34(i),1:3),YL(D34(i),1)],'-g');
end
axis equal

% % 绘制误差必须坐标一致啊
% Adata = [xy(:,1),xy(:,2),A];
% Adata = sortrows(Adata);
% Adata = sortrows(Adata,2);
% comsoldata = sortrows(comsoldata);
% comsoldata = sortrows(comsoldata,2);
% error = (comsoldata(:,3)-Adata(:,3))./comsoldata(:,3);
% Ferror = scatteredInterpolant(xy(:,1),xy(:,2),error);
% tx = 0:0.001/4:0.025;
% ty = -0.025:0.001/4:0.025;
% [qx,qy] = meshgrid(tx,ty);
% qz = Ferror(qx,qy);
% % mesh(qx,qy,qz);
% % surf(qx,qy,qz);
% figure
% contourf(qx,qy,qz,20);colorbar
% %surf(qx,qy,qz);
% title('error','FontName','Times New Roman');
% hold on
% for i=1:length(D34)
%     fill(xy(NL(D34(i),:),1),xy(NL(D34(i),:),2),'w','EdgeAlpha',0,'FaceAlpha',0.3);
%     hold on
% %     plot([XL(D34(i),1:3),XL(D34(i),1)],[YL(D34(i),1:3),YL(D34(i),1)],'-g');
% end
% axis equal

% % 绘制矩阵的右侧项
% fp = fopen('bn.txt','r');
% fgets(fp);
% Num = fscanf(fp,'%d %d',[1 2]);
% A = fscanf(fp,'%lf',[Num(1) 1]);
% fclose(fp);
% 
% F = scatteredInterpolant(xy(:,1),xy(:,2),A(1:num_nodes));
% % figure
% % F = scatteredInterpolant(X,Y,A(1:num_nodes));
% tx = 0:0.001/4:0.06;
% ty = -0.06:0.001/4:0.06;
% [qx,qy] = meshgrid(tx,ty);
% qz = F(qx,qy);
% % mesh(qx,qy,qz);
% % surf(qx,qy,qz);
% figure
% contourf(qx,qy,qz,20)
% axis equal


