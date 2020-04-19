function rotA = tetNedelec_rot( D, X, Y, Z, u)
% /*+ Function: calculate rot u of Nedelec +*/
%   int    i ;

d = 2.0 / D ;
rotA = zeros(3,1);

for i=1:6
    rotA(0+1) = rotA(0+1) + u(i) * X(i) ;
    rotA(1+1) = rotA(1+1) + u(i) * Y(i) ;
    rotA(2+1) = rotA(2+1) + u(i) * Z(i) ;
end
rotA = rotA * d;
end