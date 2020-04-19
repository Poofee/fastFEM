function a = tet_simple_a( x, y, z)
% /*+ Function: calculate determinant `a'
%                       of conventional piecewise linear tetrahedral +*/
%   /*+     | xj yj zj |
%       a = | xk yk zk |
%           | xm ym zm | +*/
%   int i;

tet_m = [0, 1, 2, 3 ;
    1, 3, 2, 0 ;
    2, 3, 0, 1 ;
    3, 1, 0, 2 ] ;
tet_m = tet_m + 1;
a = zeros(4,1);
for i=1:4
    ii = tet_m(i,1+1) ;
    jj = tet_m(i,2+1);
    kk = tet_m(i,3+1) ;
    a(i) = x(ii)*(y(jj)*z(kk) - y(kk)*z(jj))...
        +x(jj)*(y(kk)*z(ii) - y(ii)*z(kk))...
        +x(kk)*(y(ii)*z(jj) - y(jj)*z(ii)) ;
end
end