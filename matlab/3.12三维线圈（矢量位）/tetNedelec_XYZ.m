function [X, Y, Z ] = tetNedelec_XYZ( l, si, xi, yi, zi)
% /*+ Function: calculate X = l*si*xi,
%                           Y = l*si*yi,
%                           Z = l*si*zi of Nedelec +*/
X = zeros(6,1);
Y = zeros(6,1);
Z = zeros(6,1);
for i=1:6
    d = l(i) * si(i) ;
    X(i) = d * xi(i) ;
    Y(i) = d * yi(i) ;
    Z(i) = d * zi(i) ;
end
end