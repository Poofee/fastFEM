function [bdEdgeCoil] = findDomainBoundary(NL,coildomain)
% 20200327 by Poofee
% ����ÿһ������ı߽�
% Ӧ��ֻ�������޿׵�����

eCoil = int32([NL(coildomain,[1 2]); NL(coildomain,[1 3]); NL(coildomain,[2 3]);]);
[edgeCoil,i2coil,jcoil] = unique(sort(eCoil,2),'rows','legacy');
i1Coil(jcoil(size(eCoil):-1:1)) = size(eCoil):-1:1;
i1Coil = i1Coil';
bdEdgeCoilIndex = i1Coil(i1Coil == i2coil);
bdEdgeCoil = eCoil(bdEdgeCoilIndex,:);
end