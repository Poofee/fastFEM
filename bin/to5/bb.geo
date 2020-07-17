lc = 0.2;

//p1 = newp;
//Point(p1) = {0, -4.03e-3, 0, 1.0};
//+
//SetFactory("OpenCASCADE");

cH = 2.5;
cW = 0.4;
cR = 1.05;
cdx1 = 0;
cdx2 = 0;
cdy1 = -1.8;
cdy2 = -1.8;
cdz1 = 2.1;
cdz2 = -2.1;
r1 = cR+cW/2;
r2 = cR-cW/2;

//coil+
pcenter = newp;
Point(pcenter+1) = {cdx1, cdy1-cH/2, cdz1};

Point(pcenter+2) = {cdx1-r1, cdy1-cH/2, cdz1,lc};
Point(pcenter+3) = {cdx1, cdy1-cH/2, cdz1+r1,lc};
Point(pcenter+4) = {cdx1+r1, cdy1-cH/2, cdz1,lc};
Point(pcenter+5) = {cdx1, cdy1-cH/2, cdz1-r1,lc};

Point(pcenter+6) = {cdx1-r2, cdy1-cH/2, cdz1,lc};
Point(pcenter+7) = {cdx1, cdy1-cH/2, cdz1+r2,lc};
Point(pcenter+8) = {cdx1+r2, cdy1-cH/2, cdz1,lc};
Point(pcenter+9) = {cdx1, cdy1-cH/2, cdz1-r2,lc};

c5 = newc;
Circle(c5+1) = {pcenter+2,pcenter+1,pcenter+3};
Circle(c5+2) = {pcenter+3,pcenter+1,pcenter+4};
Circle(c5+3) = {pcenter+4,pcenter+1,pcenter+5};
Circle(c5+4) = {pcenter+5,pcenter+1,pcenter+2};

Circle(c5+5) = {pcenter+6,pcenter+1,pcenter+7};
Circle(c5+6) = {pcenter+7,pcenter+1,pcenter+8};
Circle(c5+7) = {pcenter+8,pcenter+1,pcenter+9};
Circle(c5+8) = {pcenter+9,pcenter+1,pcenter+6};

l5 = newll;
Curve Loop(l5) = {c5+1,c5+2,c5+3,c5+4};
l6 = newll;
Curve Loop(l6) = {c5+5,c5+6,c5+7,c5+8};
s3 = news;
Plane Surface(s3) = {l5, l6};

Extrude {0, cH, 0} {Surface{s3}; }

//coil-
pcenter = newp;
Point(pcenter+1) = {cdx2, cdy2-cH/2, cdz2};

Point(pcenter+2) = {cdx2-r1, cdy2-cH/2, cdz2,lc};
Point(pcenter+3) = {cdx2, cdy2-cH/2, cdz2+r1,lc};
Point(pcenter+4) = {cdx2+r1, cdy2-cH/2, cdz2,lc};
Point(pcenter+5) = {cdx2, cdy2-cH/2, cdz2-r1,lc};

Point(pcenter+6) = {cdx2-r2, cdy2-cH/2, cdz2,lc};
Point(pcenter+7) = {cdx2, cdy2-cH/2, cdz2+r2,lc};
Point(pcenter+8) = {cdx2+r2, cdy2-cH/2, cdz2,lc};
Point(pcenter+9) = {cdx2, cdy2-cH/2, cdz2-r2,lc};

c5 = newc;
Circle(c5+1) = {pcenter+2,pcenter+1,pcenter+3};
Circle(c5+2) = {pcenter+3,pcenter+1,pcenter+4};
Circle(c5+3) = {pcenter+4,pcenter+1,pcenter+5};
Circle(c5+4) = {pcenter+5,pcenter+1,pcenter+2};

Circle(c5+5) = {pcenter+6,pcenter+1,pcenter+7};
Circle(c5+6) = {pcenter+7,pcenter+1,pcenter+8};
Circle(c5+7) = {pcenter+8,pcenter+1,pcenter+9};
Circle(c5+8) = {pcenter+9,pcenter+1,pcenter+6};

l5 = newll;
Curve Loop(l5) = {c5+1,c5+2,c5+3,c5+4};
l6 = newll;
Curve Loop(l6) = {c5+5,c5+6,c5+7,c5+8};
s3 = news;
Plane Surface(s3) = {l5, l6};

Extrude {0, cH, 0} {Surface{s3}; }

