function [d] = tet_simple_d(x,y,z)
% /*+ Function: calculate determinant `d'
%                       of conventional piecewise linear tetrahedral +*/
%   /*+       | xj yj 1 |
%       d = - | xk yk 1 |
%             | xm ym 1 | +*/
tet_m = [0, 1, 2, 3 ;
    1, 3, 2, 0 ;
    2, 3, 0, 1 ;
    3, 1, 0, 2 ] ;
tet_m = tet_m + 1;
d = zeros(4,1);

for i=1:4
    ii = tet_m(i,1+1) ;
    jj = tet_m(i,2+1) ;
    kk = tet_m(i,3+1) ;
    d(i) = -1.0*( x(ii)*(y(jj) - y(kk))...
        +x(jj)*(y(kk) - y(ii))...
        +x(kk)*(y(ii) - y(jj)) ) ;
end
end