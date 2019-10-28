function [curlW] = RotWBasis(x,y,z,u,v,w)
% 20190830
% by Poofee
% 计算向量形函数的旋度 \nabla \times \boldsymbol{w}_j
% 棱的定义：l1=12,l2=23,l3=31,l4=14,l5=24,l6=34
edgeMap = [1,2;2,3;3,1;1,4;2,4;3,4];
curlW = zeros(6,3);
for i=1:6
    dBasis1 = dTetraNodalBasis(edgeMap(i,1),x,y,z,u,v,w);
    dBasis2 = dTetraNodalBasis(edgeMap(i,2),x,y,z,u,v,w);
    curlW(i,:) = 2*cross(dBasis1,dBasis2);
end

end