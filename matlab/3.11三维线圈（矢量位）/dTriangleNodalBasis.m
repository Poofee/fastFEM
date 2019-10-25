function [grad] = dTriangleNodalBasis(n,x,y,u,v)
% 20190828
% by Poofee
% ���������ε�Ԫ�κ������ݶ�
% n------��n���κ���
% u------��һ������
% v------�ڶ�������
if length(x) ~= 3 || length(y) ~= 3
    disp('���󣡲��������ε�Ԫ��');
end

m1 = [1,x(1),y(1);1,x(2),y(2);1,x(3),y(3)];
area = det(m1);


if n == 1
    grad = [y(2)-y(3),x(3)-x(2)]/area;
elseif n == 2
    grad = [y(3)-y(1),x(1)-x(3)]/area;
elseif n == 3
    grad = [y(1)-y(2),x(2)-x(1)]/area;
else
    grad = [];
    disp('����ı��');
end

end