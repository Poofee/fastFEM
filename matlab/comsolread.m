clear all
close all


% 获得单元和节点关系信息，并且得到每个几何Domain与其单元的关系
fp = fopen('E:\Projects\cplusplus\model\model_fine.mphtxt', 'r');
%-------------Read the head
for i=1:1:18
    tline = fgets(fp);
end
%--------------mesh point
num_nodes = fscanf(fp, '%d # number of mesh points\n', [1,1]);
pts_ind = fscanf(fp, '%d # lowest mesh point index\n', [1,1]);
fgets(fp);

xy = fscanf(fp, '%lf %lf \n', [2,num_nodes]);
xy = xy';
X=xy(:,1);
Y=xy(:,2);
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
NL=fscanf(fp, '%d %d %d \n', [3,num_elements]);
NL = NL';
NL=NL+1;
%--------------Domain
num_domain=fscanf(fp, '%d # number of geometric entity indices\n', [1,1]);
fgets(fp);
Domain = fscanf(fp, '%d \n', [1,num_domain]);
Domain = Domain';
fclose(fp);



NL(:,4) = Domain;
NDD = (1:num_nodes)';             % NDD --- 从1到ND的节点号
NDPP = NDPP + 1;
NP = length(NDPP);  % 确定边界上的节点个数，存入 NP
% 进行边界条件的设定，要通过domain进行判断，然后赋值
NDP = [];
% for i = 1:NP
% %     plot(X(NDPP(i)),Y(NDPP(i)),'b*');hold on
%     if  abs(sqrt(X(NDPP(i))^2 + Y(NDPP(i))^2) - 0.05) < 5e-3 || abs(X(NDPP(i))) < 0.5e-4
%         NDP(i) = NDPP(i);
%         plot(X(NDPP(i)),Y(NDPP(i)),'ro');
%         hold on
%     end
%
% end
NDP(NDP == 0) = [];
NDP = unique(NDP);
NP = length(NDP);  % 确定边界上的节点个数，存入 NP


num_node_free = num_nodes - NP;     % num_node_free --- 自由节点的个数

Domain3 = find(Domain==3);
Domain4 = find(Domain==4);
D34 = [Domain3;Domain4];

%%%% 设定存储的变量

CE = zeros(3,3);      % CE --- 用来存储每个单元的系数矩阵

AREA = zeros(num_elements,1);
P = zeros(num_elements,3);
Q = zeros(num_elements,3);
DD = zeros(num_elements,3);
XL = zeros(num_elements,3);
YL = zeros(num_elements,3);

colors = ['b','g','r','c','m','y','w','k'];
for i = 1:num_elements
    for j = 1:3
        k = NL(i,j);
        XL(i,j) = X(k);
        YL(i,j) = Y(k);
    end
    %     plot([XL(i,1:3),XL(i,1)],[YL(i,1:3),YL(i,1)],'-*b');
    %     hold on
    %
    %     fill(XL(i,1:3),YL(i,1:3),colors(1+mod(i,length(colors))));
    %     hold on;
    %     text(sum(XL(i,1:3))/3,sum(YL(i,1:3))/3,int2str(i));
    
    P(i,1) = YL(i,2) - YL(i,3);
    P(i,2) = YL(i,3) - YL(i,1);
    P(i,3) = YL(i,1) - YL(i,2);
    
    Q(i,1) = XL(i,3) - XL(i,2);
    Q(i,2) = XL(i,1) - XL(i,3);
    Q(i,3) = XL(i,2) - XL(i,1);
    
    DD(i,1) = XL(i,2)*YL(i,3) - XL(i,3)*YL(i,2);
    DD(i,2) = XL(i,3)*YL(i,1) - XL(i,1)*YL(i,3);
    DD(i,3) = XL(i,1)*YL(i,2) - XL(i,2)*YL(i,1);
    
    AREA(i) = 0.5 * abs(P(i,2) * Q(i,3) - Q(i,2) * P(i,3)); % AREA --- 每个单元的面积
    
    for m = 1:3
        for n = 1:3 
            CE(m,n) = (P(i,m)*P(i,n)+Q(i,m)*Q(i,n));
            
            
            if CE(m,n) > 0 && m ~= n
                fill(XL(i,1:3),YL(i,1:3),colors(2));
                hold on;
            else
                plot([XL((i),1:3),XL((i),1)],[YL((i),1:3),YL((i),1)],'-b');
                hold on;
            end
            
        end
    end
    
end
ia = vtx+1;
for i = 1:length(ia)
    plot(X(ia(i)),Y(ia(i)),'*r');
    hold on
    %text(sum(XL((ia(i)),1:3))/3,sum(YL((ia(i)),1:3))/3,int2str(ia(i)));
end
ia = entity;
for i = 1:length(ia)
    fill(XL(ia(i),1:3),YL(ia(i),1:3),'y');
    hold on
    %text(sum(XL((ia(i)),1:3))/3,sum(YL((ia(i)),1:3))/3,int2str(ia(i)));
end

for i=1:length(entity)
   plot(X(p(entity(i),1:2)),Y(p(entity(i),1:2)),'-*k') ;
   hold on
end

axis equal