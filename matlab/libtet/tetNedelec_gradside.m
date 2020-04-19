function gs = tetNedelec_gradside( D, b, c, d, Fx, Fy, Fz)
% /*+ Function: calculate (grad u, v(side)) of Nedelec +*/
%   int    i, j ;

dd = 1.0 / (6.0 * D) ;

gs = zeros(6,4);
for i=1:6
    for j=1:4
        gs(i,j) = dd * (Fx(i)*b(j) + Fy(i)*c(j) + Fz(i)*d(j)) ;
    end
end
end