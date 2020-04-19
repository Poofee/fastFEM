function [Fx, Fy, Fz] = tetNedelec_Fxyz( l, si, g, xi, yi, zi, xy, yz, zx)
% /*+ Function: calculate Fx = l*si*(yz-zi*gy+yi*gz),
%                           Fy = l*si*(zx-xi*gz+zi*gx),
%                           Fz = l*si*(xy-yi*gx+xi*gy) of Nedelec +*/
%   int    i ;

Fx = zeros(6,1);
Fy = zeros(6,1);
Fz = zeros(6,1);
for i=1:6
    d = l(i) * si(i);
    Fx(i) = d * (yz(i) - zi(i)*g(1+1) + yi(i)*g(2+1)) ;
    Fy(i) = d * (zx(i) - xi(i)*g(2+1) + zi(i)*g(0+1)) ;
    Fz(i) = d * (xy(i) - yi(i)*g(0+1) + xi(i)*g(1+1)) ;
end
end