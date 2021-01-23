function [edge,elem2edgeSign,elem2edge] = tet_ele_dof(tetmesh)
% 20200419 by Poofee
% 计算四面体单元的自由度，按照棱计算，
% edge-----所有的棱边
% elem2edgeSign-----棱边的方向，标准按照全局节点编号从低到高为1
% elem2edge------每个单元的自由度
%% 计算棱数目
NT = size(tetmesh,1);% 单元数目
% Nn = size(mesh.POS,1);% 节点数目
% 生成单元中所有的棱，顺序编号按照...
totalEdge = int32([tetmesh(:,[1 2]); tetmesh(:,[1 3]); tetmesh(:,[1 4]); ...
                   tetmesh(:,[2 3]); tetmesh(:,[2 4]); tetmesh(:,[3 4])]);
% 对边的节点编号进行处理，从低到高
sortedTotalEdge = sort(totalEdge,2);
% edge，全局的棱，大小NTx2，棱的起始点和结束点
% j，保存sortedTotalEdge中的数据在edge中的编号
[edge,i2,j] = unique(sortedTotalEdge,'rows','legacy');
Ne = size(edge,1);% 棱的数目
% elem2edge，每个单元的6条棱的全局编号，顺序为[1,2][1,3][1,4][2,3][2,4][3,4]
elem2edge = uint32(reshape(j,NT,6));
direction = ones(6*NT,1);
idx = (totalEdge(:,1)>totalEdge(:,2));
direction(idx) = -1;
% elem2edgeSign，不是编号从低到高的棱
elem2edgeSign = reshape(direction,NT,6);
end