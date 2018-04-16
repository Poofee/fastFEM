%-------------------------------------------------------------------
%                     plot things
%-------------------------------------------------------------------
fname = 'mesh.mphtxt';
[num_nodes,nodes,number_elements,nodes_ele,domain,vtx,vtx_entity,edge,tri] = read3Dmesh(fname);
%------------plot the mesh----------
X = nodes(:,1);
Y = nodes(:,2);
Z = nodes(:,3);
plot3(X,Y,Z,'*');
hold on
% plot3(X(vtx),Y(vtx),Z(vtx),'*r');
hold on
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
% for i = 1:length(nodes_ele)
%     line(X(nodes_ele(i,:)),Y(nodes_ele(i,:)),Z(nodes_ele(i,:)));
% end
%-------------plot a single tet
% n = 12;
% plot3(X(nodes_ele(n,1)),Y(nodes_ele(n,1)),Z(nodes_ele(n,1)),'*');
% hold on
% plot3(X(nodes_ele(n,[2])),Y(nodes_ele(n,[2])),Z(nodes_ele(n,[2])),'o');
% plot3(X(nodes_ele(n,[3])),Y(nodes_ele(n,[3])),Z(nodes_ele(n,[3])),'x');
% plot3(X(nodes_ele(n,[4])),Y(nodes_ele(n,[4])),Z(nodes_ele(n,[4])),'d');
% 
% plot3(X(nodes_ele(n,[1,2])),Y(nodes_ele(n,[1,2])),Z(nodes_ele(n,[1,2])));
% hold on
% plot3(X(nodes_ele(n,[2,3])),Y(nodes_ele(n,[2,3])),Z(nodes_ele(n,[2,3])));
% plot3(X(nodes_ele(n,[3,4])),Y(nodes_ele(n,[3,4])),Z(nodes_ele(n,[3,4])));
% plot3(X(nodes_ele(n,[1,3])),Y(nodes_ele(n,[1,3])),Z(nodes_ele(n,[1,3])));
% plot3(X(nodes_ele(n,[1,4])),Y(nodes_ele(n,[1,4])),Z(nodes_ele(n,[1,4])));
% plot3(X(nodes_ele(n,[2,4])),Y(nodes_ele(n,[2,4])),Z(nodes_ele(n,[2,4])));
%------------check the direction of the tet
%in comsol, it follows the right hand rule
% for i = 1:length(nodes_ele)
%     K = [X(nodes_ele(i,1)),Y(nodes_ele(i,1)),Z(nodes_ele(i,1))];
%     M = [X(nodes_ele(i,2)),Y(nodes_ele(i,2)),Z(nodes_ele(i,2))];
%     N = [X(nodes_ele(i,3)),Y(nodes_ele(i,3)),Z(nodes_ele(i,3))];
%     L = [X(nodes_ele(i,4)),Y(nodes_ele(i,4)),Z(nodes_ele(i,4))];
%     KM = M - K;
%     MN = N - M;
%     KL = L - K;
%     c = cross(KM,MN);
%     d(i) = dot(c,KL);
% end
axis equal