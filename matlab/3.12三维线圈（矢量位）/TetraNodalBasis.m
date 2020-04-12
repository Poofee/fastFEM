function [basis] = TetraNodalBasis(n,x,y,z,u,v,w)
% 20190828
% by Poofee
% һ���Ľڵ������嵥Ԫ�κ���
% n-------��n���κ���,n=1,2,3,4
% u-------��1������
% v-------��2������
% w-------��3������
% basis------����õ����κ�����ֵ
if length(x) ~= 4 || length(y) ~= 4 || length(z) ~= 4
    disp('���󣡲����Ľڵ������嵥Ԫ��');
end

m1 = [reshape(x,[4,1]),reshape(y,[4,1]),reshape(z,[4,1])];
V = det([1;1;1]*m1(4,:)-m1(1:3,:));
m1(n,:) = [u,v,w];
Vn = det([1;1;1]*m1(4,:)-m1(1:3,:));
basis = Vn/V;

end