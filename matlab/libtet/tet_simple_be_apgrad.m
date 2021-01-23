function [ag] = tet_simple_be_apgrad(D, b, c, d, u)
dd = 1.0 / 24.0 ;

ag = zeros(3,1);
for i=1:4
    ag(i) = 0.0 ;
    for j=1:4
        ag(i) = ag(i) + (b(i)*u(j,1) + c(i)*u(j,2) + d(i)*u(j,3)) * dd ;
    end
end
end