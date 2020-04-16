function tetg = tet_Center( x, y, z)
% /*+ Function: calculate center of gravity of tetrahedral +*/
tetg = zeros(3,1);
tetg(1) = ( x(4) + x(1) + x(2) + x(3) ) * 0.25 ;
tetg(2) = ( y(4) + y(1) + y(2) + y(3) ) * 0.25 ;
tetg(3) = ( z(4) + z(1) + z(2) + z(3) ) * 0.25 ;
end