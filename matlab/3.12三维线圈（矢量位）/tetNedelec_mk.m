function [xi,yi,zi] = tetNedelec_mk(x,y,z)
% Function: calculate xm-xk, ym-yk, zm-zk of Nedelec
tet_n = [ 0 , 1 , 2 , 3  ;
    0 , 2 , 3 , 1 ;
    0 , 3 , 1 , 2 ;
    1 , 2 , 0 , 3 ;
    1 , 3 , 2 , 0 ;
    2 , 3 , 0 , 1 ];
tet_n = tet_n + 1;
xi = zeros(6,1);
yi = zeros(6,1);
zi = zeros(6,1);
for i=1:6
    ii = tet_n(i,2+1) ;
    jj = tet_n(i,3+1) ;
    xi(i) = x(jj) - x(ii) ;
    yi(i) = y(jj) - y(ii) ;
    zi(i) = z(jj) - z(ii) ;
end
end