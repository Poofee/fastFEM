function [grad] = dTriangleNodalBasis(n,x,y,u,v)
% 20190828
% by Poofee
% 计算三角形单元形函数的梯度
% n------第n个形函数
% u------第一个坐标
% v------第二个坐标
if length(x) ~= 3 || length(y) ~= 3
    disp('错误！不是三角形单元。');
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
    disp('错误的编号');
end

end