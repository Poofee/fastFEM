function[gg] = tet_simple_gradgrad(D, b, c, d)
% Function: calculate (grad u, grad v) of tetrahedral
dd = 1.0 / (6.0 * D) ;
gg = zeros(4,4);

for i=1:4
    for j=1:4
        gg(i,j) = (b(i)*b(j) + c(i)*c(j) + d(i)*d(j)) * dd ;
    end
end

end