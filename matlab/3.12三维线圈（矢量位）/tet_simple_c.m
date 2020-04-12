function [c] = tet_simple_c(x,y,z)
% /*+ Function: calculate determinant `c'
%                       of conventional piecewise linear tetrahedral +*/
%   /*+       | xj 1 zj |
%       c = - | xk 1 zk |
%             | xm 1 zm | +*/
tet_m = [0, 1, 2, 3 ;
    1, 3, 2, 0 ;
    2, 3, 0, 1 ;
    3, 1, 0, 2 ] ;
tet_m = tet_m + 1;
c = zeros(4,1);

for i=1:4
    ii = tet_m(i,1+1) ;
    jj = tet_m(i,2+1) ;
    kk = tet_m(i,3+1) ;
    c(i) = -1.0*( z(ii)*(x(jj) - x(kk))...
        +z(jj)*(x(kk) - x(ii))...
        +z(kk)*(x(ii) - x(jj)) ) ;
end
end