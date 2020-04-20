function [grad] = dTetraNodalBasis(x,y,z,u,v,w)
% 20190828
% by Poofee
% 计算四节点四面体单元形函数的梯度
% n------第n个形函数,n=1,2,3,4
% u------第一个坐标，对于一阶单元来说可能没用，因为是常量
% v------第二个坐标
% w------第三个坐标
if length(x) ~= 4 || length(y) ~= 4 || length(z) ~= 4
    disp('错误！不是四节点四面体单元。');
end
grad = zeros(4,3);
m1 = [reshape(x,[4,1]),reshape(y,[4,1]),reshape(z,[4,1])];
V6 = det([1;1;1]*m1(4,:)-m1(1:3,:));
if V6 < 0
    fprintf('V6小于零\n');
end
% if n == 1
%     VdBdx = [-1       ,0        ,0        ;
%              x(4)-x(2),y(4)-y(2),z(4)-z(2);
%              x(4)-x(3),y(4)-y(3),z(4)-z(3)];
%     VdBdy = [0        ,-1       ,0        ;
%              x(4)-x(2),y(4)-y(2),z(4)-z(2);
%              x(4)-x(3),y(4)-y(3),z(4)-z(3)];
%     VdBdz = [0        ,0        ,-1       ;
%              x(4)-x(2),y(4)-y(2),z(4)-z(2);
%              x(4)-x(3),y(4)-y(3),z(4)-z(3)];
%     grad = [det(VdBdx),det(VdBdy),det(VdBdz)]/V6;
    grad(1,:) = cross([x(3)-x(2),y(3)-y(2),z(3)-z(2)],[x(2)-x(4),y(2)-y(4),z(2)-z(4)])/V6;
% elseif n == 2
%     VdBdx = [x(4)-x(1),y(4)-y(1),z(4)-z(1);
%              -1       ,0         ,0        ;
%              x(4)-x(3),y(4)-y(3),z(4)-z(3)];
%     VdBdy = [x(4)-x(1),y(4)-y(1),z(4)-z(1);
%              0        ,-1       ,0        ;
%              x(4)-x(3),y(4)-y(3),z(4)-z(3)];
%     VdBdz = [x(4)-x(1),y(4)-y(1),z(4)-z(1);
%              0        ,0        ,-1       ;
%              x(4)-x(3),y(4)-y(3),z(4)-z(3)];
%     grad = [det(VdBdx),det(VdBdy),det(VdBdz)]/V6;
    grad(2,:) = cross([x(3)-x(1),y(3)-y(1),z(3)-z(1)],[x(4)-x(1),y(4)-y(1),z(4)-z(1)])/V6;
% elseif n == 3
%     VdBdx = [x(4)-x(1),y(4)-y(1),z(4)-z(1);
%              x(4)-x(2),y(4)-y(2),z(4)-z(2);
%              -1       ,0        ,0       ];
%     VdBdy = [x(4)-x(1),y(4)-y(1),z(4)-z(1);
%              x(4)-x(2),y(4)-y(2),z(4)-z(2);
%              0        ,-1       ,0       ];
%     VdBdz = [x(4)-x(1),y(4)-y(1),z(4)-z(1);
%              x(4)-x(2),y(4)-y(2),z(4)-z(2);
%              0        ,0        ,-1      ];
%     grad = [det(VdBdx),det(VdBdy),det(VdBdz)]/V6;
    grad(3,:) = cross([x(4)-x(1),y(4)-y(1),z(4)-z(1)],[x(2)-x(1),y(2)-y(1),z(2)-z(1)])/V6;
% elseif n == 4
%     grad = -(dTetraNodalBasis(1,x,y,z,0,0,0)...
%             +dTetraNodalBasis(2,x,y,z,0,0,0)...
%             +dTetraNodalBasis(3,x,y,z,0,0,0));
    grad(4,:) = cross([x(2)-x(1),y(2)-y(1),z(2)-z(1)],[x(3)-x(1),y(3)-y(1),z(3)-z(1)])/V6;
% else
%     grad = [];
%     disp('节点编号错误');
% end

end