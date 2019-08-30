function [W] = WBasis(x,y,z,u,v,w)
% 20190830
% by Poofee
% 计算向量形函数
% 棱的定义：l1=12,l2=23,l3=31,l4=14,l5=24,l6=34
edgeMap = [1,2;2,3;3,1;1,4;2,4;3,4];
W = zeros(6,1);
for i=1:6
  W(i,1) = TetraNodalBasis(edgeMap(i,1),x,y,z,u,v,w)*dTetraNodalBasis(edgeMap(i,2),x,y,z,u,v,w)...
    -TetraNodalBasis(edgeMap(i,2),x,y,z,u,v,w)*dTetraNodalBasis(edgeMap(i,1),x,y,z,u,v,w);  
end

end