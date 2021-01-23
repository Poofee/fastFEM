function [isBdEdge,isBdNode,freeEdge,freeNode] = findDomain3DBoundary(meshTETS,edge,elem2edge)
% 20200419 by Poofee
% 计算分网区域的边界
% 输入参数：
% meshTETS-----------区域内的分网单元，可以是某个domain的，也可以是所有的
% edge---------------全局的棱变量
% elem2edge-----------全局的单元内的棱变量
% 输出参数：
% isBdEdge------------棱的边界信息
% isBdNode------------边界节点

% 生成四面体的所有面
allFace = [meshTETS(:,[2 4 3]);meshTETS(:,[1 3 4]);meshTETS(:,[1 4 2]);meshTETS(:,[1 2 3])];

Nfall = length(allFace);% 面的数目
% 不重复的face
[face, i2, j] = unique(sort(allFace,2),'rows','legacy');
i1(j(Nfall:-1:1)) = Nfall:-1:1;% 只有出现一次的编号才不会变
i1 = i1';
bdFlag = zeros(Nfall,1,'uint8');
bdFaceidx = i1(i1==i2);% 边界上的face的索引
bdFlag(bdFaceidx) = 1;% Dirichlet边界条件
bdFlag = reshape(bdFlag,size(meshTETS,1),4);% 转换为对应的单元形式
% 查找边界棱和节点
isBdEdge = false(size(edge,1),1);
isBdNode = true(max(edge,[],'all'),1);% 感觉应该是找到最大的编号，前提是，编号从1开始，顺序递增
% 由于编号问题，POS变量保存的点并不全是分网节点
allNode = meshTETS(:,1:4);
allNode = allNode(:);
% uniNode = unique(allNode);
isBdNode(allNode) = false;% 找到全部分网节点
% 将边界面上的点设为边界点
% isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
% isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
% isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
% isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
% bdEdge = edge(isBdEdge,:);
% isBdNode(bdEdge) = true;% 标记边界节点
% % 绘制边界，验证
% hold on
% for i=1:size(bdEdge,1)
%     line(mesh.POS(bdEdge(i,:),1),mesh.POS(bdEdge(i,:),2),mesh.POS(bdEdge(i,:),3),'Color',[0 0 0]);
% end
% DisplayNodes(mesh.POS(isBdNode,:));
% view(45,atan(1/sqrt(2))*180/pi);
% 未知的变量
freeEdge = find(~isBdEdge);
freeNode = find(~isBdNode);
fprintf('自由棱数目为%d,自由节点数目为%d\n',size(freeEdge,1),size(freeNode,1));
end