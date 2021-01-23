function [isBdEdge,isBdNode,freeEdge,freeNode] = findDomain3DBoundary(meshTETS,edge,elem2edge)
% 20200419 by Poofee
% �����������ı߽�
% ���������
% meshTETS-----------�����ڵķ�����Ԫ��������ĳ��domain�ģ�Ҳ���������е�
% edge---------------ȫ�ֵ������
% elem2edge-----------ȫ�ֵĵ�Ԫ�ڵ������
% ���������
% isBdEdge------------��ı߽���Ϣ
% isBdNode------------�߽�ڵ�

% �����������������
allFace = [meshTETS(:,[2 4 3]);meshTETS(:,[1 3 4]);meshTETS(:,[1 4 2]);meshTETS(:,[1 2 3])];

Nfall = length(allFace);% �����Ŀ
% ���ظ���face
[face, i2, j] = unique(sort(allFace,2),'rows','legacy');
i1(j(Nfall:-1:1)) = Nfall:-1:1;% ֻ�г���һ�εı�ŲŲ����
i1 = i1';
bdFlag = zeros(Nfall,1,'uint8');
bdFaceidx = i1(i1==i2);% �߽��ϵ�face������
bdFlag(bdFaceidx) = 1;% Dirichlet�߽�����
bdFlag = reshape(bdFlag,size(meshTETS,1),4);% ת��Ϊ��Ӧ�ĵ�Ԫ��ʽ
% ���ұ߽���ͽڵ�
isBdEdge = false(size(edge,1),1);
isBdNode = true(max(edge,[],'all'),1);% �о�Ӧ�����ҵ����ı�ţ�ǰ���ǣ���Ŵ�1��ʼ��˳�����
% ���ڱ�����⣬POS��������ĵ㲢��ȫ�Ƿ����ڵ�
allNode = meshTETS(:,1:4);
allNode = allNode(:);
% uniNode = unique(allNode);
isBdNode(allNode) = false;% �ҵ�ȫ�������ڵ�
% ���߽����ϵĵ���Ϊ�߽��
% isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
% isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
% isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
% isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
% bdEdge = edge(isBdEdge,:);
% isBdNode(bdEdge) = true;% ��Ǳ߽�ڵ�
% % ���Ʊ߽磬��֤
% hold on
% for i=1:size(bdEdge,1)
%     line(mesh.POS(bdEdge(i,:),1),mesh.POS(bdEdge(i,:),2),mesh.POS(bdEdge(i,:),3),'Color',[0 0 0]);
% end
% DisplayNodes(mesh.POS(isBdNode,:));
% view(45,atan(1/sqrt(2))*180/pi);
% δ֪�ı���
freeEdge = find(~isBdEdge);
freeNode = find(~isBdNode);
fprintf('��������ĿΪ%d,���ɽڵ���ĿΪ%d\n',size(freeEdge,1),size(freeNode,1));
end