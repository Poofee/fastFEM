function [basis] = TriangleNodalBasis(n,x,y,u,v)
% 20190828
% by Poofee
% 一阶三节点三角形单元形函数
% n-------第n个形函数
% u-------第1个坐标
% v-------第2个坐标
% w-------第3个坐标
% basis------计算得到的形函数的值
if length(x) ~= 3 || length(y) ~= 3
    disp('错误！不是三角形单元。');
end

m1 = [1,x(1),y(1);1,x(2),y(2);1,x(3),y(3)];
area = det(m1);

m1(n,:) = [1,u,v];
basis = det(m1)/area;

end