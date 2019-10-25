function [basis] = NodalBasisFunctions(n,element,x,y,z,u,v,w)
% 20190828
% by Poofee
% 节点有限元形函数，参数含义：
% n-----第n个形函数
% element------单元类型
% u-------计算点的第一个坐标
% v-------计算点的第二个坐标
% w-------计算点的第三个坐标
% basis------计算得到的形函数的值

if element == 303% 三节点三角形单元
    basis = TriangleNodalBasis(n,x,y,u,v);
elseif element == 404% 四节点四面体单元
    basis = TetraNodalBasis(n,x,y,z,u,v,w);
else
    basis = [];
    disp('错误，没有找到对应的单元类型！');
end

end
