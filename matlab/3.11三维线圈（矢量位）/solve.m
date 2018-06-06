% 2018-04-14
% by Poofee
% This program is a test for 3D magnetistaic FEM code.
% The model in this code is simple coil.
% �ο��������̲ġ��������̵�ų���ֵ��������89-108ҳ
% ��ȡ�����ļ�0
close all
clear all
fname = 'mesh.mphtxt';
[num_nodes,nodes,number_elements,nodes_ele,domain] = read3Dmesh(fname);

X = nodes(:,1);
Y = nodes(:,2);
Z = nodes(:,3);
%-----------����һЩ��������
%-----------��Ԫ���
volume = zeros(length(nodes_ele),1);
for i = 1:length(nodes_ele)
   v_tmp = [ones(4,1),X(nodes_ele(i,:)),Y(nodes_ele(i,:)),Z(nodes_ele(i,:))]; 
   volume(i) = det(v_tmp) / 6;
end
%-----------���� pK pM pN pL
p = zeros(length(nodes_ele),4);
for i = 1:length(nodes_ele)
   pK = [X(nodes_ele(i,[2 3 4])),Y(nodes_ele(i,[2 3 4])),Z(nodes_ele(i,[2 3 4]))]; 
   pM = -[X(nodes_ele(i,[1 3 4])),Y(nodes_ele(i,[1 3 4])),Z(nodes_ele(i,[1 3 4]))];
   pN = [X(nodes_ele(i,[1 2 4])),Y(nodes_ele(i,[1 2 4])),Z(nodes_ele(i,[1 2 4]))];
   pL = -[X(nodes_ele(i,[1 2 3])),Y(nodes_ele(i,[1 2 3])),Z(nodes_ele(i,[1 2 3]))];
   p(i,:) = [det(pK),det(pM),det(pN),det(pL)];
end
%-----------���� qK qM qN qL
q = zeros(length(nodes_ele),4);
for i = 1:length(nodes_ele)
   qK = -[ones(3,1),Y(nodes_ele(i,[2 3 4])),Z(nodes_ele(i,[2 3 4]))]; 
   qM = [ones(3,1),Y(nodes_ele(i,[1 3 4])),Z(nodes_ele(i,[1 3 4]))];
   qN = -[ones(3,1),Y(nodes_ele(i,[1 2 4])),Z(nodes_ele(i,[1 2 4]))];
   qL = [ones(3,1),Y(nodes_ele(i,[1 2 3])),Z(nodes_ele(i,[1 2 3]))];
   q(i,:) = [det(qK),det(qM),det(qN),det(qL)];
end
%-----------���� rK rM rN rL
r = zeros(length(nodes_ele),4);
for i = 1:length(nodes_ele)
   rK = -[X(nodes_ele(i,[2 3 4])),ones(3,1),Z(nodes_ele(i,[2 3 4]))]; 
   rM = [X(nodes_ele(i,[1 3 4])),ones(3,1),Z(nodes_ele(i,[1 3 4]))];
   rN = -[X(nodes_ele(i,[1 2 4])),ones(3,1),Z(nodes_ele(i,[1 2 4]))];
   rL = [X(nodes_ele(i,[1 2 3])),ones(3,1),Z(nodes_ele(i,[1 2 3]))];
   r(i,:) = [det(rK),det(rM),det(rN),det(rL)];
end
%-----------���� sK sM sN sL
s = zeros(length(nodes_ele),4);
for i = 1:length(nodes_ele)
   sK = -[X(nodes_ele(i,[2 3 4])),Y(nodes_ele(i,[2 3 4])),ones(3,1)]; 
   sM = [X(nodes_ele(i,[1 3 4])),Y(nodes_ele(i,[1 3 4])),ones(3,1)];
   sN = -[X(nodes_ele(i,[1 2 4])),Y(nodes_ele(i,[1 2 4])),ones(3,1)];
   sL = [X(nodes_ele(i,[1 2 3])),Y(nodes_ele(i,[1 2 3])),ones(3,1)];
   s(i,:) = [det(sK),det(sM),det(sN),det(sL)];
end

% ��������Ԫ�����
% ��Ҫע����ǣ���������ʸ�����������ά����3n

% ��ȡ�߽�

S = zeros(3*num_nodes,3*num_nodes);
eS = zeros(3,3);

% ����Ԫװ��
% for i = 1:number_elements
%     for j = 1:4
%     	IR = nodes_ele(i,j);
%         IR = find(arrangenode==IR);
%         for k = 1:4
%         	IN = nodes_ele(i,k);
%             IN = find(arrangenode==IN);
%             
%             eS(1,1) = (r(i,j)*r(i,k)+s(i,j)*s(i,k)*mu/36/volume(i);
%             eS(1,2) = -(r(i,j)*q(i,k))*mu/36/volume(i);
%             eS(1,3) = -(s(i,j)*q(i,k))*mu/36/volume(i);
% 
%             eS(2,1) = -(q(i,j)*r(i,k))*mu/36/volume(i);
%             eS(2,2) = (q(i,j)*q(i,k)+s(i,j)*s(i,k)*mu/36/volume(i);
%             eS(2,3) = -(s(i,j)*r(i,k))*mu/36/volume(i);
% 
%             eS(3,1) = -(q(i,j)*s(i,k))*mu/36/volume(i);
%             eS(3,2) = -(r(i,j)*s(i,k))*mu/36/volume(i);
%             eS(3,3) = (q(i,j)*q(i,k)+r(i,j)*r(i,k)*mu/36/volume(i);
% 
%             S(IR,IN) = S(IR,IN) + eS(j,k);
%         end
% 
%     end
% end

