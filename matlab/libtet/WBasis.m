function [W] = WBasis(x,y,z,u,v,w)
% 20190830
% by Poofee
% ���������κ���
% ��Ķ��壺l1=12,l2=23,l3=31,l4=14,l5=24,l6=34
edgeMap = [1,2;1,3;1,4;2,3;2,4;3,4];
W = zeros(6,3);
gradN = dTetraNodalBasis(x,y,z,u,v,w);
for i=1:6
%     L = norm([x(edgeMap(i,1))-x(edgeMap(i,2)),y(edgeMap(i,1))-y(edgeMap(i,2)),z(edgeMap(i,1))-z(edgeMap(i,2))]);
    W(i,:) = TetraNodalBasis(edgeMap(i,1),x,y,z,u,v,w)*gradN(edgeMap(i,2),:)...
        -TetraNodalBasis(edgeMap(i,2),x,y,z,u,v,w)*gradN(edgeMap(i,1),:);
    %   W(i,:) = W(i,:)/L;
end

end