function [] = setBoundary3D(meshTet,boudValue)
fprintf('开始设置边界条件...\n');
% 一般求解的时候，不需要特别的设置边界条件，因为有一个默认的几何边界，那么为了
% 设置边界，就需要找出这个边界。以前的思路是拿几何尺寸来判断，缺点就是你必须知道
% 边界的形状和尺寸，非常不方便。但是，再想一下，如果一个单元是边界，二维的时候，
% 这个棱边界只存在于一个分网单元中，那么就可以查找出现一次的边，那就是边界了。
% 三维的思路也是一样的，只不过要查找的是面单元。

% 设置Neumann边界条件

% 设置Dirichlet边界条件
u = zeros(Ne,1);
if ~isempty(edge(isBdEdge,:))
    u(isBdEdge) = 0;% 固定边界条件
    %     f = f - (A - M)*u;
    f(isBdEdge) = u(isBdEdge);
end
g0 = -M'*u;
end