lc = 4;
nn = 40; // mesh subdivisions per turn

DefineConstant
[
 ht = {46, Name "Parameters/Tube height"},
 rt1 = {16, Name "Parameters/Tube internal radius"},
 rt2 = {24, Name "Parameters/Tube external radius"},
 lb = {150, Name "Parameters/Infinite box width"}
];

// tube
p = newp;
Point(p) = {0, 0, -ht/2, lc};
Point(p+1) = {rt1, 0, -ht/2, lc};
Point(p+2) = {0, rt1, -ht/2, lc};
Point(p+3) = {-rt1, 0, -ht/2, lc};
Point(p+4) = {0, -rt1, -ht/2, lc};
Point(p+5) = {rt2, 0, -ht/2, lc};
Point(p+6) = {0, rt2, -ht/2, lc};
Point(p+7) = {-rt2, 0, -ht/2, lc};
Point(p+8) = {0, -rt2, -ht/2, lc};
c = newl;
Circle(c) = {p+1, p, p+2};
Circle(c+1) = {p+2, p, p+3};
Circle(c+2) = {p+3, p, p+4};
Circle(c+3) = {p+4, p, p+1};
Circle(c+4) = {p+5, p, p+6};
Circle(c+5) = {p+6, p, p+7};
Circle(c+6) = {p+7, p, p+8};
Circle(c+7) = {p+8, p, p+5};
ll = newll;
Curve Loop(ll) = {c+4, c+5, c+6, c+7, -c, -(c+1), -(c+2), -(c+3)};
s = news;
Plane Surface(s) = {ll};
tmp[] = Extrude {0,0,ht}{ Surface{s}; /*Layers{nn};*/ };
vol_tube = tmp[1];

// box
p = newp;
Point(p) = {-lb/2,-lb/2,-lb/2};
Point(p+1) = {lb/2,-lb/2,-lb/2};
Point(p+2) = {lb/2,lb/2,-lb/2};
Point(p+3) = {-lb/2,lb/2,-lb/2};
Point(p+4) = {-lb/2,-lb/2,lb/2};
Point(p+5) = {lb/2,-lb/2,lb/2};
Point(p+6) = {lb/2,lb/2,lb/2};
Point(p+7) = {-lb/2,lb/2,lb/2};

l = newl;
Line(l) = {p,p+1};
Line(l+1) = {p+1,p+2};
Line(l+2) = {p+2,p+3};
Line(l+3) = {p+3,p};
Line(l+4) = {p+4,p+5};
Line(l+5) = {p+5,p+6};
Line(l+6) = {p+6,p+7};
Line(l+7) = {p+7,p+4};
Line(l+8) = {p, p+4};
Line(l+9) = {p+1, p+5};
Line(l+10) = {p+2, p+6};
Line(l+11) = {p+3, p+7};

ll = newll;
Curve Loop(ll+1) = {l+8, -(l+7), -(l+11), l+3};
Curve Loop(ll+2) = {l+9, l+5, -(l+10), -(l+1)};
Curve Loop(ll+3) = {l,l+1,l+2,l+3};
Curve Loop(ll+4) = {l+4,l+5,l+6,l+7};
Curve Loop(ll+5) = {l+2, l+11, -(l+6), -(l+10)};
Curve Loop(ll+6) = {l, l+9, -(l+4), -(l+8)};

s = news;
Plane Surface(s) = {ll+1};
Plane Surface(s+1) = {ll+2};
Plane Surface(s+2) = {ll+3};
Plane Surface(s+3) = {ll+4};
Plane Surface(s+4) = {ll+5};
Plane Surface(s+5) = {ll+6};

sl = newsl;
Surface Loop(sl) = {s:s+5};
Surface Loop(sl+1) = CombinedBoundary{ Volume{vol_tube[]}; };

v = newv;
Volume(v) = {sl, sl+1};