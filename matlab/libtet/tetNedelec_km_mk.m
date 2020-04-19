function [xy, yz, zx ] = tetNedelec_km_mk( x, y, z)
% /*+ Function: calculate xkym-xmyk, ykzm-ymzk, zkxm-zmxk of Nedelec +*/
%   int    i ;
tet_n = [ 0 , 1 , 2 , 3  ;
    0 , 2 , 3 , 1 ;
    0 , 3 , 1 , 2 ;
    1 , 2 , 0 , 3 ;
    1 , 3 , 2 , 0 ;
    2 , 3 , 0 , 1 ];
tet_n = tet_n + 1;
xy = zeros(6,1);
yz = zeros(6,1);
zx = zeros(6,1);
for i=1:6
    ii = tet_n(i,2+1);
    jj = tet_n(i,3+1);
    xy(i) = x(ii)*y(jj) - x(jj)*y(ii) ;
    yz(i) = y(ii)*z(jj) - y(jj)*z(ii) ;
    zx(i) = z(ii)*x(jj) - z(jj)*x(ii) ;
end
end