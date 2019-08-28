function [basis] = NodalBasisFunctions(n,element,x,y,z,u,v,w)
% 20190828
% by Poofee
% �ڵ�����Ԫ�κ������������壺
% n-----��n���κ���
% element------��Ԫ����
% u-------�����ĵ�һ������
% v-------�����ĵڶ�������
% w-------�����ĵ���������
% basis------����õ����κ�����ֵ

if element == 303% ���ڵ������ε�Ԫ
    basis = TriangleNodalBasis(n,x,y,u,v);
elseif element == 404% �Ľڵ������嵥Ԫ
    basis = TetraNodalBasis(n,x,y,z,u,v,w);
else
    basis = [];
    disp('����û���ҵ���Ӧ�ĵ�Ԫ���ͣ�');
end

end
