function [basis] = TriangleNodalBasis(n,x,y,u,v)
% 20190828
% by Poofee
% һ�����ڵ������ε�Ԫ�κ���
% n-------��n���κ���
% u-------��1������
% v-------��2������
% w-------��3������
% basis------����õ����κ�����ֵ
if length(x) ~= 3 || length(y) ~= 3
    disp('���󣡲��������ε�Ԫ��');
end

m1 = [1,x(1),y(1);1,x(2),y(2);1,x(3),y(3)];
area = det(m1);

m1(n,:) = [1,u,v];
basis = det(m1)/area;

end