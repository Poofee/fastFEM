function [be_gs] = tetNedelec_be_gradside(u,D,b,c,d,Fx,Fy,Fz)
%   /*+ Function: calculate (grad u, v(side)) of Nedelec +*/

be_gs = zeros(6,1);

gs = tetNedelec_gradside( D, b, c, d, Fx, Fy, Fz) ;
for i=1:6
    be_gs(i) = 0.0 ;
    for j=1:4
        be_gs(i) = be_gs(i) + u(j) * gs(i,j) ;
    end
end
end