function [b] = tet_simple_b( x, y, z)
% /*+ Function: calculate determinant `b'
%                       of conventional piecewise linear tetrahedral +*/
%   /*+       | 1 yj zj |
%       b = - | 1 yk zk |
%             | 1 ym zm | +*/
tet_m = [0, 1, 2, 3 ;
    1, 3, 2, 0 ;
    2, 3, 0, 1 ;
    3, 1, 0, 2 ] ;
tet_m = tet_m + 1;
b = zeros(4,1);
for i=1:4
    ii = tet_m(i,1+1) ;
    jj = tet_m(i,2+1) ;
    kk = tet_m(i,3+1) ;
    b(i) = -1.0*( y(ii)*(z(jj) - z(kk))...
        +y(jj)*(z(kk) - z(ii))...
        +y(kk)*(z(ii) - z(jj)) ) ;
end
end