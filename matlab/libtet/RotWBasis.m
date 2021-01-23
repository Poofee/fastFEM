function [curlW] = RotWBasis(x,y,z,u,v,w,locEdge)
% 20190830
% by Poofee
% 计算向量形函数的旋度 \nabla \times \boldsymbol{w}_j
% 棱的定义：l1=12,l2=23,l3=31,l4=14,l5=24,l6=34
% edgeMap = [1,2;2,3;3,1;1,4;2,4;3,4];
% edgeMap = [1,2;1,3;1,4;2,3;2,4;3,4];
edgeMap = locEdge;%[1 2;2 3;3 1;1 4;2 4;3 4;];
curlW = zeros(6,3);
gradN = dTetraNodalBasis(x,y,z,u,v,w);
for i=1:6
    dBasis1 = gradN(edgeMap(i,1),:);
    dBasis2 = gradN(edgeMap(i,2),:);
    curlW(i,:) = 2*cross(dBasis1,dBasis2);
end

end