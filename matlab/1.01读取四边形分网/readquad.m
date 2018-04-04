clear all
close all

dims = 2;
nnodes = 4;
% 获得单元和节点关系信息，并且得到每个几何Domain与其单元的关系
fp = fopen('..\model\reg1.mphtxt', 'r');
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
NL(:,5) = Domain;

XL = zeros(num_elements,4);
YL = zeros(num_elements,4);
colors = ['b','g','r','c','m','y','w','k'];
for i = 1:num_elements
    for j = 1:4
        k = NL(i,j);
        XL(i,j) = xy(k,1);
        YL(i,j) = xy(k,2);        
        text(XL(i,j),YL(i,j),num2str(j));
        hold on
    end
    plot(XL(i,:),YL(i,:),'-*');        
    hold off
        %pause(1)
end


axis equal