function [basis] = TetraNodalBasis(n,x,y,z,u,v,w)
% 20190828
% by Poofee
% 一阶四节点四面体单元形函数
% n-------第n个形函数,n=1,2,3,4
% u-------第1个坐标
% v-------第2个坐标
% w-------第3个坐标
% basis------计算得到的形函数的值
if length(x) ~= 4 || length(y) ~= 4 || length(z) ~= 4
    disp('错误！不是四节点四面体单元。');
end

m1 = [reshape(x,[4,1]),reshape(y,[4,1]),reshape(z,[4,1])];
V = det([1;1;1]*m1(4,:)-m1(1:3,:));
m1(n,:) = [u,v,w];
Vn = det([1;1;1]*m1(4,:)-m1(1:3,:));
basis = Vn/V;

end