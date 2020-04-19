function [msh] = ReadCOMSOL3D(fileName)
% 20191026 by Poofee
% 读取三维COMSOL分网文件
% 2016-09-25
% by Poofee
% 2016-10-20
% by Poofee
msh.Types = { ...
    { 2, 1, 'LINES', 'nbLines'}, ... % 1
    { 3,  2, 'TRIANGLES', 'nbTriangles'}, ...
    { 4,  2, 'QUADS', 'nbQuads'}, ...  
    { 4,  3, 'TETS', 'nbTets'}, ...
    { 8,  3, 'HEXAS', 'nbHexas'}, ... %5
    { 6,  3, 'PRISMS', 'nbPrisms'}, ...
    { 5,  3, 'PYRAMIDS', 'nbPyramids'}, ...
    { 3,  1, 'LINES3', 'nbLines3'}, ...
    { 6,  2, 'TRIANGLES6', 'nbTriangles6'}, ...
    { 9,  2, 'QUADS9', 'nbQuads9'}, ... % 10
    { 10,  3, 'TETS10', 'nbTets10'}, ...
    { 27,  3, 'HEXAS27', 'nbHexas27'}, ...
    { 18,  3, 'PRISMS18', 'nbPrisms18'}, ...
    { 14,  3, 'PYRAMIDS14', 'nbPyramids14'}, ...
    { 1,  0, 'POINTS', 'nbPoints'}, ... % 15
    { 8,  3, 'QUADS8', 'nbQuads8'}, ...
    { 20,  3, 'HEXAS20', 'nbHexas20'}, ...
    { 15,  3, 'PRISMS15', 'nbPrisms15'}, ...
    { 13,  3, 'PYRAMIDS13', 'nbPyramids13'}, ...
};
%-------------------------------------------------------------------
%                     deal with the file
%-------------------------------------------------------------------
fp = fopen(fileName, 'r');
if fp < 0
   disp(['错误：打开文件',fileName,'失败']); 
end
%-------------Read the head
for i=1:1:17
    tline = fgets(fp);
end
%--------------mesh point
num_nodes = fscanf(fp, '%d # number of mesh points\n', [1,1]);
msh.nbNod = num_nodes;
pts_ind = fscanf(fp, '%d # lowest mesh point index\n', [1,1]);
fgets(fp);

xyz = fscanf(fp, '%lf %lf \n', [3,num_nodes]);
xyz = xyz';
% X = xyz(:,1);
% Y = xyz(:,2);
% Z = xyz(:,3); 
msh.POS = xyz;
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
% plot3(X,Y,Z,'*');
%--------------vertex
for i=1:7
    fgets(fp);
end
num_vtx_ns = fscanf(fp, '%d # number of nodes per element\n', [1,1]);
num_vtx_ele = fscanf(fp, '%d # number of elements\n', [1,1]);
fgets(fp);
vtx = fscanf(fp, '%d \n', [1,num_vtx_ele]);
vtx = vtx' + 1;
% plot3(X(vtx),Y(vtx),Z(vtx),'or');hold on
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
%%--------------plot the line
% for i = 1:length(edge)
%    plot3(X(edge(i,:)),Y(edge(i,:)),Z(edge(i,:)),'-*'); hold on
%    text(X(edge(i,1)),Y(edge(i,1)),Z(edge(i,1)),num2str(edge_entity(i)-1));
%    text(X(edge(i,2)),Y(edge(i,2)),Z(edge(i,2)),num2str(edge_entity(i)-1));
% end
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
msh.nbTets = number_elements;
fgets(fp);
fgets(fp);
nodes_ele = fscanf(fp, '%d %d %d %d\n', [4,number_elements]);
nodes_ele = nodes_ele' + 1;%index starts from zero
msh.TETS = nodes_ele;
%-------------domain
domain_num = fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
domain = fscanf(fp,'%d \n',[1,domain_num]);
domain_num = domain_num';
msh.DOMAIN = domain';
fclose(fp);

%-------------------------------------------------------------------
%                     plot things
%-------------------------------------------------------------------
%------------plot the mesh----------
% plot3(X,Y,Z,'*');
% hold on
% plot3(X(vtx),Y(vtx),Z(vtx),'*r');
% hold on
% for i = 1:length(vtx_entity)
%     n = vtx_entity(i);
%     plot3(X(nodes_ele(n,[1,2])),Y(nodes_ele(n,[1,2])),Z(nodes_ele(n,[1,2])));
%     hold on
%     plot3(X(nodes_ele(n,[2,3])),Y(nodes_ele(n,[2,3])),Z(nodes_ele(n,[2,3])));
%     plot3(X(nodes_ele(n,[3,4])),Y(nodes_ele(n,[3,4])),Z(nodes_ele(n,[3,4])));
%     plot3(X(nodes_ele(n,[1,3])),Y(nodes_ele(n,[1,3])),Z(nodes_ele(n,[1,3])));
%     plot3(X(nodes_ele(n,[1,4])),Y(nodes_ele(n,[1,4])),Z(nodes_ele(n,[1,4])));
%     plot3(X(nodes_ele(n,[2,4])),Y(nodes_ele(n,[2,4])),Z(nodes_ele(n,[2,4])));
% end
% for i = 1:length(nodes_element)
%     line(X(nodes_element(i,:)),Y(nodes_element(i,:)),Z(nodes_element(i,:)));
% end
% %-------------plot a single tet
% n = 12;
% plot3(X(nodes_element(n,1)),Y(nodes_element(n,1)),Z(nodes_element(n,1)),'*');
% hold on
% plot3(X(nodes_element(n,[2])),Y(nodes_element(n,[2])),Z(nodes_element(n,[2])),'o');
% plot3(X(nodes_element(n,[3])),Y(nodes_element(n,[3])),Z(nodes_element(n,[3])),'x');
% plot3(X(nodes_element(n,[4])),Y(nodes_element(n,[4])),Z(nodes_element(n,[4])),'d');
% 
% plot3(X(nodes_element(n,[1,2])),Y(nodes_element(n,[1,2])),Z(nodes_element(n,[1,2])));
% hold on
% plot3(X(nodes_element(n,[2,3])),Y(nodes_element(n,[2,3])),Z(nodes_element(n,[2,3])));
% plot3(X(nodes_element(n,[3,4])),Y(nodes_element(n,[3,4])),Z(nodes_element(n,[3,4])));
% plot3(X(nodes_element(n,[1,3])),Y(nodes_element(n,[1,3])),Z(nodes_element(n,[1,3])));
% plot3(X(nodes_element(n,[1,4])),Y(nodes_element(n,[1,4])),Z(nodes_element(n,[1,4])));
% plot3(X(nodes_element(n,[2,4])),Y(nodes_element(n,[2,4])),Z(nodes_element(n,[2,4])));
% %------------check the direction of the tet
% %in comsol, it follows the right hand rule
% for i = 1:length(nodes_element)
%     K = [X(nodes_element(i,1)),Y(nodes_element(i,1)),Z(nodes_element(i,1))];
%     M = [X(nodes_element(i,2)),Y(nodes_element(i,2)),Z(nodes_element(i,2))];
%     N = [X(nodes_element(i,3)),Y(nodes_element(i,3)),Z(nodes_element(i,3))];
%     L = [X(nodes_element(i,4)),Y(nodes_element(i,4)),Z(nodes_element(i,4))];
%     KM = M - K;
%     MN = N - M;
%     KL = L - K;
%     c = cross(KM,MN);
%     d(i) = dot(c,KL);
% end
end
























