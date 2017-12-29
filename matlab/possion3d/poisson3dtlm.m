% 2016-10-21
% by Poofee
% tlm version
% 
% use this program solve a parallal cap
% clear all;
% close all;
%------------------some define
epsilon = 8.85e-12;
%-------------------------------------------------------------------
%                     deal with the file
%-------------------------------------------------------------------
fp = fopen('xihua.mphtxt', 'r');
%-------------Read the head
for i=1:1:18
    tline = fgets(fp);
end
%--------------mesh point
num_nodes = fscanf(fp, '%d # number of mesh points\n', [1,1]);
pts_ind = fscanf(fp, '%d # lowest mesh point index\n', [1,1]);
fgets(fp);

xyz = fscanf(fp, '%lf %lf \n', [3,num_nodes]);
xyz = xyz';
X = xyz(:,1);
Y = xyz(:,2);
Z = xyz(:,3); 
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
fgets(fp);
fgets(fp);
nodes_ele = fscanf(fp, '%d %d %d %d\n', [4,number_elements]);
nodes_ele = nodes_ele' + 1;%index starts from zero
%-------------domain
domain_num = fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
domain = fscanf(fp,'%d \n',[1,domain_num]);
domain_num = domain_num';
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

%-------------------------------------------------------------------
%                     electrostatic 3d
%K = 1 M = 2 N = 3 L = 4
%-------------------------------------------------------------------
%-----------volume of single element
volume = zeros(length(nodes_ele),1);
for i = 1:length(nodes_ele)
   v_tmp = [ones(4,1),X(nodes_ele(i,:)),Y(nodes_ele(i,:)),Z(nodes_ele(i,:))]; 
   volume(i) = det(v_tmp) / 6;
end
%-----------calculate pK pM pN pL
p = zeros(length(nodes_ele),4);
for i = 1:length(nodes_ele)
   pK = [X(nodes_ele(i,[2 3 4])),Y(nodes_ele(i,[2 3 4])),Z(nodes_ele(i,[2 3 4]))]; 
   pM = -[X(nodes_ele(i,[1 3 4])),Y(nodes_ele(i,[1 3 4])),Z(nodes_ele(i,[1 3 4]))];
   pN = [X(nodes_ele(i,[1 2 4])),Y(nodes_ele(i,[1 2 4])),Z(nodes_ele(i,[1 2 4]))];
   pL = -[X(nodes_ele(i,[1 2 3])),Y(nodes_ele(i,[1 2 3])),Z(nodes_ele(i,[1 2 3]))];
   p(i,:) = [det(pK),det(pM),det(pN),det(pL)];
end
%-----------calculate qK qM qN qL
q = zeros(length(nodes_ele),4);
for i = 1:length(nodes_ele)
   qK = -[ones(3,1),Y(nodes_ele(i,[2 3 4])),Z(nodes_ele(i,[2 3 4]))]; 
   qM = [ones(3,1),Y(nodes_ele(i,[1 3 4])),Z(nodes_ele(i,[1 3 4]))];
   qN = -[ones(3,1),Y(nodes_ele(i,[1 2 4])),Z(nodes_ele(i,[1 2 4]))];
   qL = [ones(3,1),Y(nodes_ele(i,[1 2 3])),Z(nodes_ele(i,[1 2 3]))];
   q(i,:) = [det(qK),det(qM),det(qN),det(qL)];
end
%-----------calculate rK rM rN rL
r = zeros(length(nodes_ele),4);
for i = 1:length(nodes_ele)
   rK = -[X(nodes_ele(i,[2 3 4])),ones(3,1),Z(nodes_ele(i,[2 3 4]))]; 
   rM = [X(nodes_ele(i,[1 3 4])),ones(3,1),Z(nodes_ele(i,[1 3 4]))];
   rN = -[X(nodes_ele(i,[1 2 4])),ones(3,1),Z(nodes_ele(i,[1 2 4]))];
   rL = [X(nodes_ele(i,[1 2 3])),ones(3,1),Z(nodes_ele(i,[1 2 3]))];
   r(i,:) = [det(rK),det(rM),det(rN),det(rL)];
end
%-----------calculate sK sM sN sL
s = zeros(length(nodes_ele),4);
for i = 1:length(nodes_ele)
   sK = -[X(nodes_ele(i,[2 3 4])),Y(nodes_ele(i,[2 3 4])),ones(3,1)]; 
   sM = [X(nodes_ele(i,[1 3 4])),Y(nodes_ele(i,[1 3 4])),ones(3,1)];
   sN = -[X(nodes_ele(i,[1 2 4])),Y(nodes_ele(i,[1 2 4])),ones(3,1)];
   sL = [X(nodes_ele(i,[1 2 3])),Y(nodes_ele(i,[1 2 3])),ones(3,1)];
   s(i,:) = [det(sK),det(sM),det(sN),det(sL)];
end
%--------------check the point if on the boundary
Ra = 0.02; 
Rb = 0.03;

for i = 1:num_nodes
   radius = X.^2 + Y.^2;
   radius = sqrt(radius);
   boundarya = find(radius < (0.02 + 1e-5));
   boundaryb = find(radius > (0.03 - 1e-5));
   unknown = find(abs(radius - 0.025) < 0.004);
end
arrangenode = [unknown;boundarya;boundaryb];
%--------------assemble the global matrix
S = zeros(num_nodes,num_nodes);
eS = zeros(4,4);
for i = 1:number_elements
    for j = 1:4
        for k = 1:4
            eS(j,k) = (q(i,j)*q(i,k)+r(i,j)*r(i,k)+s(i,j)*s(i,k))*epsilon/36/volume(i);
        end
    end
    for j = 1:4
        IR = nodes_ele(i,j);
        IR = find(arrangenode==IR);
        for k = 1:4
            IN = nodes_ele(i,k);
            IN = find(arrangenode==IN);
            S(IR,IN) = S(IR,IN) + eS(j,k);
        end
    end
end

%--------------rearrange the global matrix no need

%--------------block matrix
%       [C11 C12] [Vx] = [Fc]
%       [C12 C22] [Vc] = [Fx]
%       C11Vx + C12Vc = Fc
%       Vx = C11 \ (Fc - C12Vc)  ;
C11 = S(1:length(unknown),1:length(unknown));
C12 = S(1:length(unknown),(length(unknown)+1):length(S));
V = zeros(num_nodes,1);
Vx = V(1:length(unknown));
Vc = V((length(unknown)+1):length(S));
Vc(1:length(boundarya)) = 10;
Vc((length(boundarya)+1):length(Vc)) = 0;
Vx = C11 \(-C12 * Vc);      
plot(100*radius(unknown),Vx,'*')

phi = 10 * log(Rb./radius(unknown))/log(Rb/Ra);
% rr = [2.0833 2.1666 2.25 2.333 2.4166 2.5 2.5833 2.666 2.75 2.83 2.9166];
% phi = [8.9932097 8.0259115 7.0951129 6.19818 5.3327146 4.4966029 3.6879095 2.904884 2.1459646 1.4097028 0.694785];
hold on
plot(radius(unknown),phi,'ro')

error = abs(phi - Vx)./phi;
























