function [edge,elem2edgeSign,elem2edge] = tet_ele_dof(tetmesh)
% 20200419 by Poofee
% ���������嵥Ԫ�����ɶȣ���������㣬
% edge-----���е����
% elem2edgeSign-----��ߵķ��򣬱�׼����ȫ�ֽڵ��Ŵӵ͵���Ϊ1
% elem2edge------ÿ����Ԫ�����ɶ�
%% ��������Ŀ
NT = size(tetmesh,1);% ��Ԫ��Ŀ
% Nn = size(mesh.POS,1);% �ڵ���Ŀ
% ���ɵ�Ԫ�����е��⣬˳���Ű���...
totalEdge = int32([tetmesh(:,[1 2]); tetmesh(:,[1 3]); tetmesh(:,[1 4]); ...
                   tetmesh(:,[2 3]); tetmesh(:,[2 4]); tetmesh(:,[3 4])]);
% �ԱߵĽڵ��Ž��д����ӵ͵���
sortedTotalEdge = sort(totalEdge,2);
% edge��ȫ�ֵ��⣬��СNTx2�������ʼ��ͽ�����
% j������sortedTotalEdge�е�������edge�еı��
[edge,i2,j] = unique(sortedTotalEdge,'rows','legacy');
Ne = size(edge,1);% �����Ŀ
% elem2edge��ÿ����Ԫ��6�����ȫ�ֱ�ţ�˳��Ϊ[1,2][1,3][1,4][2,3][2,4][3,4]
elem2edge = uint32(reshape(j,NT,6));
direction = ones(6*NT,1);
idx = (totalEdge(:,1)>totalEdge(:,2));
direction(idx) = -1;
% elem2edgeSign�����Ǳ�Ŵӵ͵��ߵ���
elem2edgeSign = reshape(direction,NT,6);
end