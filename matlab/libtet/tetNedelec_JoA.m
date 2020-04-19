function [be_JoA] = tetNedelec_JoA( Joe, si, l, b, c, d)
%     /*+ Function: calculate (Jo(apex*3), A*) of Nedelec +*/
tet_n = [ 0 , 1 , 2 , 3  ;
    0 , 2 , 3 , 1 ;
    0 , 3 , 1 , 2 ;
    1 , 2 , 0 , 3 ;
    1 , 3 , 2 , 0 ;
    2 , 3 , 0 , 1 ];
tet_n = tet_n + 1;

be_JoA = zeros(6,1);
PN = zeros(6,12);
d0 = [2.0,  1.0,  1.0,  1.0;
    2.0,  1.0,  1.0,  1.0;
    2.0,  1.0,  1.0,  1.0;
    1.0,  2.0,  1.0,  1.0;
    1.0,  2.0,  1.0,  1.0;
    1.0,  1.0,  2.0,  1.0];
d1 = [ -1.0, -2.0, -1.0, -1.0;
    -1.0, -1.0, -2.0, -1.0;
    -1.0, -1.0, -1.0, -2.0;
    -1.0, -1.0, -2.0, -1.0;
    -1.0, -1.0, -1.0, -2.0;
    -1.0, -1.0, -1.0, -2.0;];

for i=1:6
    dd = si(i) * l(i) / 120.0 ;
    for j=1:12
        PN(i,j) = dd ;
    end
end

for i=1:6
    for j=1:4
        PN(i,j) = PN(i,j) * (d0(i,j)*b(tet_n(i,1+1)) + d1(i,j)*b(tet_n(i,0+1)));
    end
end
for i=1:6
    k=4;
    for j=1:4
        k = k+1;
        PN(i,k) = PN(i,k)*(d0(i,j)*c(tet_n(i,1+1)) + d1(i,j)*c(tet_n(i,0+1)));
    end
end
for i=1:6
    k = 8;
    for j=1:4
        k = k + 1;
        PN(i,k) = PN(i,k)*(d0(i,j)*d(tet_n(i,1+1)) + d1(i,j)*d(tet_n(i,0+1)));
    end
end

for i=1:6
    be_JoA(i) = 0.0 ;
    for j=1:4
        be_JoA(i) = be_JoA(i)+PN(i,j) * Joe(j,0+1) ;
    end
end
for i=1:6
    k=4;
    for j=1:4
        k=k+1;
        be_JoA(i) = be_JoA(i)+PN(i,k) * Joe(j,1+1) ;
    end
end
for i=1:6
    k=8;
    for j=1:4
        k=k+1;
        be_JoA(i) = be_JoA(i)+PN(i,k) * Joe(j,2+1) ;
    end
end
end