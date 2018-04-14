function [num_nodes,nodes,number_elements,nodes_ele,domain] = read3Dmesh(fname)
%READ3DMESH 此处显示有关此函数的摘要
%   此处显示详细说明
fp = fopen(fname 'r');
%-------------读取文件头
for i=1:1:18
    tline = fgets(fp);
end
%--------------读取分网节点坐标
num_nodes = fscanf(fp, '%d # number of mesh points\n', [1,1]);
pts_ind = fscanf(fp, '%d # lowest mesh point index\n', [1,1]);
fgets(fp);

xyz = fscanf(fp, '%lf %lf \n', [3,num_nodes]);
xyz = xyz';
X = xyz(:,1);
Y = xyz(:,2);
Z = xyz(:,3); 
nodes = [X;Y;Z];
%--------------vertex
for i=1:7
    fgets(fp);
end
num_vtx_ns = fscanf(fp, '%d # number of nodes per element\n', [1,1]);
num_vtx_ele = fscanf(fp, '%d # number of elements\n', [1,1]);
fgets(fp);
vtx = fscanf(fp, '%d \n', [1,num_vtx_ele]);
vtx = vtx' + 1;
%--------------vertex entity
num_vtx_entity = fscanf(fp, '%d # number of geometric entity indices\n', [1,1]);
fgets(fp);
vtx_entity = fscanf(fp, '%d \n', [1,num_vtx_entity]);
vtx_entity = vtx_entity' + 1;
%--------------edge
for i=1:5
    fgets(fp);
end
num_edge_ns = fscanf(fp, '%d # number of nodes per element\n', [1,1]);
num_edge = fscanf(fp, '%d # number of elements\n', [1,1]);
fgets(fp);
edge = fscanf(fp, '%d %d\n', [2,num_edge]);
edge = edge' + 1;
%--------------edge entity
num_edge_entity=fscanf(fp, '%d # number of geometric entity indices\n', [1,1]);
fgets(fp);
edge_entity = fscanf(fp, '%d \n', [1,num_edge_entity]);
edge_entity = edge_entity' + 1;
%--------------tri
for i=1:5
    fgets(fp);
end
ns_per_ele=fscanf(fp, '%d # number of nodes per element\n', [1,1]);
num_tri=fscanf(fp, '%d # number of elements\n', [1,1]);
fgets(fp);
tri = fscanf(fp, '%d %d %d \n', [3,num_tri]);
tri = tri' + 1;
%--------------tri entity
num_tri_entity = fscanf(fp, '%d # number of geometric entity indices\n', [1,1]);
fgets(fp);
tri_entity = fscanf(fp, '%d \n', [1,num_tri_entity]);
tri_entity = tri_entity' + 1;
%-------------tet
for i=1:5
    fgets(fp);
end
num_per_ele = fscanf(fp,'%d # number of nodes per element',[1,1]);
number_elements = fscanf(fp,'%d # number of elements',[1,1]);
fgets(fp);
fgets(fp);
nodes_ele = fscanf(fp, '%d %d %d %d\n', [4,number_elements]);
nodes_ele = nodes_ele' + 1;%index starts from zero
%-------------domain
domain_num = fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
domain = fscanf(fp,'%d \n',[1,domain_num]);
domain = domain';
fclose(fp);

end

