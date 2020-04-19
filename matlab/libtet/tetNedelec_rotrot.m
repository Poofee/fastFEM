function [rr] = tetNedelec_rotrot( D, X, Y, Z)
% /*+ Function: calculate (rot u, rot v) of Nedelec +*/

d = 2.0 / (3.0 * D) ;

rr = zeros(6,6);
for i=1:6
    for j=1:6
        rr(i,j) = d * (X(i)*X(j) + Y(i)*Y(j) + Z(i)*Z(j)) ;
    end
end
end