//SetFactory("OpenCASCADE");

halfWidth = 8.05;
//down part
//left
p1 = newp;
Point(p1+1) = {0, -halfWidth, -34};
Point(p1+2) = {0, -halfWidth, -12.15};
Point(p1+3) = {0.98, -halfWidth, -12.15};
Point(p1+4) = {4.825, -halfWidth, 2.2};
Point(p1+5) = {5.95, -halfWidth, 2.2};
Point(p1+6) = {5.95, -halfWidth, -28};
Point(p1+7) = {17.1, -halfWidth, -28};
Point(p1+8) = {17.1, -halfWidth, 22.25};
Point(p1+9) = {17.686, -halfWidth, 22.25};
Point(p1+10) = {23, -halfWidth, 7.65};
Point(p1+11) = {25.85, -halfWidth, 7.65};
Point(p1+12) = {25.85, -halfWidth, 5.475};
Point(p1+13) = {23, -halfWidth, 5.475};
Point(p1+14) = {23, -halfWidth, 3.3};
Point(p1+15) = {25.85, -halfWidth, 3.3};
Point(p1+16) = {25.85, -halfWidth, -2.557};
Point(p1+17) = {23.25, -halfWidth, -9.7};
Point(p1+18) = {23.25, -halfWidth, -34};

l1=newl;
//Line (l1+1) = {p1+1, p1+2} ;
Line (l1+2) = {p1+2, p1+3} ;
Line (l1+3) = {p1+3, p1+4} ;
Line (l1+4) = {p1+4, p1+5} ;
Line (l1+5) = {p1+5, p1+6} ;
Line (l1+6) = {p1+6, p1+7} ;
Line (l1+7) = {p1+7, p1+8} ;
Line (l1+8) = {p1+8, p1+9} ;
Line (l1+9) = {p1+9, p1+10} ;
Line (l1+10) = {p1+10, p1+11} ;
Line (l1+11) = {p1+11, p1+12} ;
Line (l1+12) = {p1+12, p1+13} ;
Line (l1+13) = {p1+13, p1+14} ;
Line (l1+14) = {p1+14, p1+15} ;
Line (l1+15) = {p1+15, p1+16} ;
Line (l1+16) = {p1+16, p1+17} ;
Line (l1+17) = {p1+17, p1+18} ;
Line (l1+18) = {p1+18, p1+1} ;

//duichen
p11 = newp;
//Point(p11+1) = {0, -halfWidth, -34};
//Point(p11+2) = {0, -halfWidth, -12.15};
Point(p11+3) = {-0.98, -halfWidth, -12.15};
Point(p11+4) = {-4.825, -halfWidth, 2.2};
Point(p11+5) = {-5.95, -halfWidth, 2.2};
Point(p11+6) = {-5.95, -halfWidth, -28};
Point(p11+7) = {-17.1, -halfWidth, -28};
Point(p11+8) = {-17.1, -halfWidth, 22.25};
Point(p11+9) = {-17.686, -halfWidth, 22.25};
Point(p11+10) = {-23, -halfWidth, 7.65};
Point(p11+11) = {-25.85, -halfWidth, 7.65};
Point(p11+12) = {-25.85, -halfWidth, 5.475};
Point(p11+13) = {-23, -halfWidth, 5.475};
Point(p11+14) = {-23, -halfWidth, 3.3};
Point(p11+15) = {-25.85, -halfWidth, 3.3};
Point(p11+16) = {-25.85, -halfWidth, -2.557};
Point(p11+17) = {-23.25, -halfWidth, -9.7};
Point(p11+18) = {-23.25, -halfWidth, -34};

l11=newl;
//Line (l11+1) = {p11+1, p11+2} ;
//Line (l11+2) = {p11+2, p11+3} ;
Line (l11+3) = {p11+3, p11+4} ;
Line (l11+4) = {p11+4, p11+5} ;
Line (l11+5) = {p11+5, p11+6} ;
Line (l11+6) = {p11+6, p11+7} ;
Line (l11+7) = {p11+7, p11+8} ;
Line (l11+8) = {p11+8, p11+9} ;
Line (l11+9) = {p11+9, p11+10} ;
Line (l11+10) = {p11+10, p11+11} ;
Line (l11+11) = {p11+11, p11+12} ;
Line (l11+12) = {p11+12, p11+13} ;
Line (l11+13) = {p11+13, p11+14} ;
Line (l11+14) = {p11+14, p11+15} ;
Line (l11+15) = {p11+15, p11+16} ;
Line (l11+16) = {p11+16, p11+17} ;
Line (l11+17) = {p11+17, p11+18} ;
//Line (l11+18) = {p11+18, p11+1} ;

Line (l11+19) = {p11+18, p1+1} ;
Line (l11+20) = {p11+3, p1+2} ;

//right
p1r = newp;
Point(p1r+1) = {0, halfWidth, -34};
Point(p1r+2) = {0, halfWidth, -12.15};
Point(p1r+3) = {0.98, halfWidth, -12.15};
Point(p1r+4) = {4.825, halfWidth, 2.2};
Point(p1r+5) = {5.95, halfWidth, 2.2};
Point(p1r+6) = {5.95, halfWidth, -28};
Point(p1r+7) = {17.1, halfWidth, -28};
Point(p1r+8) = {17.1, halfWidth, 22.25};
Point(p1r+9) = {17.686, halfWidth, 22.25};
Point(p1r+10) = {23, halfWidth, 7.65};
Point(p1r+11) = {25.85, halfWidth, 7.65};
Point(p1r+12) = {25.85, halfWidth, 5.475};
Point(p1r+13) = {23, halfWidth, 5.475};
Point(p1r+14) = {23, halfWidth, 3.3};
Point(p1r+15) = {25.85, halfWidth, 3.3};
Point(p1r+16) = {25.85, halfWidth, -2.557};
Point(p1r+17) = {23.25, halfWidth, -9.7};
Point(p1r+18) = {23.25, halfWidth, -34};

l1r=newl;
//Line (l1r+1) = {p1r+1, p1r+2} ;
Line (l1r+2) = {p1r+2, p1r+3} ;
Line (l1r+3) = {p1r+3, p1r+4} ;
Line (l1r+4) = {p1r+4, p1r+5} ;
Line (l1r+5) = {p1r+5, p1r+6} ;
Line (l1r+6) = {p1r+6, p1r+7} ;
Line (l1r+7) = {p1r+7, p1r+8} ;
Line (l1r+8) = {p1r+8, p1r+9} ;
Line (l1r+9) = {p1r+9, p1r+10} ;
Line (l1r+10) = {p1r+10, p1r+11} ;
Line (l1r+11) = {p1r+11, p1r+12} ;
Line (l1r+12) = {p1r+12, p1r+13} ;
Line (l1r+13) = {p1r+13, p1r+14} ;
Line (l1r+14) = {p1r+14, p1r+15} ;
Line (l1r+15) = {p1r+15, p1r+16} ;
Line (l1r+16) = {p1r+16, p1r+17} ;
Line (l1r+17) = {p1r+17, p1r+18} ;
Line (l1r+18) = {p1r+18, p1r+1} ;
//right
p11r = newp;
//Point(p11r+1) = {0, halfWidth, -34};
//Point(p11r+2) = {0, halfWidth, -12.15};
Point(p11r+3) = {-0.98, halfWidth, -12.15};
Point(p11r+4) = {-4.825, halfWidth, 2.2};
Point(p11r+5) = {-5.95, halfWidth, 2.2};
Point(p11r+6) = {-5.95, halfWidth, -28};
Point(p11r+7) = {-17.1, halfWidth, -28};
Point(p11r+8) = {-17.1, halfWidth, 22.25};
Point(p11r+9) = {-17.686, halfWidth, 22.25};
Point(p11r+10) = {-23, halfWidth, 7.65};
Point(p11r+11) = {-25.85, halfWidth, 7.65};
Point(p11r+12) = {-25.85, halfWidth, 5.475};
Point(p11r+13) = {-23, halfWidth, 5.475};
Point(p11r+14) = {-23, halfWidth, 3.3};
Point(p11r+15) = {-25.85, halfWidth, 3.3};
Point(p11r+16) = {-25.85, halfWidth, -2.557};
Point(p11r+17) = {-23.25, halfWidth, -9.7};
Point(p11r+18) = {-23.25, halfWidth, -34};

l11r=newl;
//Line (l11r+1) = {p11r+1, p11r+2} ;
//Line (l11r+2) = {p11r+2, p11r+3} ;
Line (l11r+3) = {p11r+3, p11r+4} ;
Line (l11r+4) = {p11r+4, p11r+5} ;
Line (l11r+5) = {p11r+5, p11r+6} ;
Line (l11r+6) = {p11r+6, p11r+7} ;
Line (l11r+7) = {p11r+7, p11r+8} ;
Line (l11r+8) = {p11r+8, p11r+9} ;
Line (l11r+9) = {p11r+9, p11r+10} ;
Line (l11r+10) = {p11r+10, p11r+11} ;
Line (l11r+11) = {p11r+11, p11r+12} ;
Line (l11r+12) = {p11r+12, p11r+13} ;
Line (l11r+13) = {p11r+13, p11r+14} ;
Line (l11r+14) = {p11r+14, p11r+15} ;
Line (l11r+15) = {p11r+15, p11r+16} ;
Line (l11r+16) = {p11r+16, p11r+17} ;
Line (l11r+17) = {p11r+17, p11r+18} ;
//Line (l11+18) = {p11+18, p11+1} ;

Line (l11r+19) = {p11r+18, p1r+1} ;
Line (l11r+20) = {p11r+3, p1r+2} ;

//l2 = newll;
//Curve Loop(l2) = {l1+1, l1+2, l1+3, l1+4, l1+5, l1+6, l1+7, l1+8, l1+9, l1+10, l1+11, l1+12, l1+13, l1+14, l1+15, l1+16, l1+17, l1+18};

//s1 = news;
//Plane Surface(s1) = {l2};

lm1 = newl;
Line (lm1+1) = {p1+1,p1r+1} ;
Line (lm1+2) = {p1+2,p1r+2} ;
Line (lm1+3) = {p1+3,p1r+3} ;
Line (lm1+4) = {p1+4,p1r+4} ;
Line (lm1+5) = {p1+5,p1r+5} ;
Line (lm1+6) = {p1+6,p1r+6} ;
Line (lm1+7) = {p1+7,p1r+7} ;
Line (lm1+8) = {p1+8,p1r+8} ;
Line (lm1+9) = {p1+9,p1r+9} ;
Line (lm1+10) = {p1+10,p1r+10} ;
Line (lm1+11) = {p1+11,p1r+11} ;
Line (lm1+12) = {p1+12,p1r+12} ;
Line (lm1+13) = {p1+13,p1r+13} ;
Line (lm1+14) = {p1+14,p1r+14} ;
Line (lm1+15) = {p1+15,p1r+15} ;
Line (lm1+16) = {p1+16,p1r+16} ;
Line (lm1+17) = {p1+17,p1r+17} ;
Line (lm1+18) = {p1+18,p1r+18} ;

lm2 = newl;
Line (lm2+3) = {p11+3,p11r+3} ;
Line (lm2+4) = {p11+4,p11r+4} ;
Line (lm2+5) = {p11+5,p11r+5} ;
Line (lm2+6) = {p11+6,p11r+6} ;
Line (lm2+7) = {p11+7,p11r+7} ;
Line (lm2+8) = {p11+8,p11r+8} ;
Line (lm2+9) = {p11+9,p11r+9} ;
Line (lm2+10) = {p11+10,p11r+10} ;
Line (lm2+11) = {p11+11,p11r+11} ;
Line (lm2+12) = {p11+12,p11r+12} ;
Line (lm2+13) = {p11+13,p11r+13} ;
Line (lm2+14) = {p11+14,p11r+14} ;
Line (lm2+15) = {p11+15,p11r+15} ;
Line (lm2+16) = {p11+16,p11r+16} ;
Line (lm2+17) = {p11+17,p11r+17} ;
Line (lm2+18) = {p11+18,p11r+18} ;

//up part
//right
p2 = newp;
Point(p2+1) = {0, halfWidth, -11.15};
Point(p2+2) = {0.595, halfWidth, -11.15};
Point(p2+3) = {4.44, halfWidth, 3.2};
Point(p2+4) = {5.5, halfWidth, 3.2};
Point(p2+5) = {5.5, halfWidth, 30};
Point(p2+6) = {18, halfWidth, 30};
Point(p2+7) = {18, halfWidth, 23.625};
Point(p2+8) = {23.45, halfWidth, 8.65};
Point(p2+9) = {25.35, halfWidth, 8.65};
Point(p2+10) = {25.35, halfWidth, 25};
Point(p2+11) = {23.2, halfWidth, 25};
Point(p2+12) = {23.2, halfWidth, 35};
Point(p2+13) = {0, halfWidth, 35};

l3 = newl;
Line (l3+1) = {p2+1, p2+2} ;
Line (l3+2) = {p2+2, p2+3} ;
Line (l3+3) = {p2+3, p2+4} ;
Line (l3+4) = {p2+4, p2+5} ;
Line (l3+5) = {p2+5, p2+6} ;
Line (l3+6) = {p2+6, p2+7} ;
Line (l3+7) = {p2+7, p2+8} ;
Line (l3+8) = {p2+8, p2+9} ;
Line (l3+9) = {p2+9, p2+10} ;
Line (l3+10) = {p2+10, p2+11} ;
Line (l3+11) = {p2+11, p2+12} ;
Line (l3+12) = {p2+12, p2+13} ;
//Line (l3+13) = {p2+13, p2+1} ;

//duichen
p22 = newp;
//Point(p22+1) = {0, halfWidth, -11.15};
Point(p22+2) = {-0.595, halfWidth, -11.15};
Point(p22+3) = {-4.44, halfWidth, 3.2};
Point(p22+4) = {-5.5, halfWidth, 3.2};
Point(p22+5) = {-5.5, halfWidth, 30};
Point(p22+6) = {-18, halfWidth, 30};
Point(p22+7) = {-18, halfWidth, 23.625};
Point(p22+8) = {-23.45, halfWidth, 8.65};
Point(p22+9) = {-25.35, halfWidth, 8.65};
Point(p22+10) = {-25.35, halfWidth, 25};
Point(p22+11) = {-23.2, halfWidth, 25};
Point(p22+12) = {-23.2, halfWidth, 35};
//Point(p22+13) = {-0, halfWidth, 35};

l33 = newl;
//Line (l33+1) = {p22+1, p22+2} ;
Line (l33+2) = {p22+2, p22+3} ;
Line (l33+3) = {p22+3, p22+4} ;
Line (l33+4) = {p22+4, p22+5} ;
Line (l33+5) = {p22+5, p22+6} ;
Line (l33+6) = {p22+6, p22+7} ;
Line (l33+7) = {p22+7, p22+8} ;
Line (l33+8) = {p22+8, p22+9} ;
Line (l33+9) = {p22+9, p22+10} ;
Line (l33+10) = {p22+10, p22+11} ;
Line (l33+11) = {p22+11, p22+12} ;
//Line (l33+12) = {p22+12, p22+13} ;
//Line (l33+13) = {p22+13, p22+1} ;

Line (l33+14) = {p22+2, p2+1} ;
Line (l33+15) = {p22+12, p2+13} ;

p2r = newp;
Point(p2r+1) = {0, -halfWidth, -11.15};
Point(p2r+2) = {0.595, -halfWidth, -11.15};
Point(p2r+3) = {4.44, -halfWidth, 3.2};
Point(p2r+4) = {5.5, -halfWidth, 3.2};
Point(p2r+5) = {5.5, -halfWidth, 30};
Point(p2r+6) = {18, -halfWidth, 30};
Point(p2r+7) = {18, -halfWidth, 23.625};
Point(p2r+8) = {23.45, -halfWidth, 8.65};
Point(p2r+9) = {25.35, -halfWidth, 8.65};
Point(p2r+10) = {25.35, -halfWidth, 25};
Point(p2r+11) = {23.2, -halfWidth, 25};
Point(p2r+12) = {23.2, -halfWidth, 35};
Point(p2r+13) = {0, -halfWidth, 35};

l3r = newl;
Line (l3r+1) = {p2r+1, p2r+2} ;
Line (l3r+2) = {p2r+2, p2r+3} ;
Line (l3r+3) = {p2r+3, p2r+4} ;
Line (l3r+4) = {p2r+4, p2r+5} ;
Line (l3r+5) = {p2r+5, p2r+6} ;
Line (l3r+6) = {p2r+6, p2r+7} ;
Line (l3r+7) = {p2r+7, p2r+8} ;
Line (l3r+8) = {p2r+8, p2r+9} ;
Line (l3r+9) = {p2r+9, p2r+10} ;
Line (l3r+10) = {p2r+10, p2r+11} ;
Line (l3r+11) = {p2r+11, p2r+12} ;
Line (l3r+12) = {p2r+12, p2r+13} ;
//Line (l3r+13r) = {p2r+13r, p2r+1} ;

//right
p22r = newp;
//Point(p22r+1) = {0, -halfWidth, -11.15};
Point(p22r+2) = {-0.595, -halfWidth, -11.15};
Point(p22r+3) = {-4.44, -halfWidth, 3.2};
Point(p22r+4) = {-5.5, -halfWidth, 3.2};
Point(p22r+5) = {-5.5, -halfWidth, 30};
Point(p22r+6) = {-18, -halfWidth, 30};
Point(p22r+7) = {-18, -halfWidth, 23.625};
Point(p22r+8) = {-23.45, -halfWidth, 8.65};
Point(p22r+9) = {-25.35, -halfWidth, 8.65};
Point(p22r+10) = {-25.35, -halfWidth, 25};
Point(p22r+11) = {-23.2, -halfWidth, 25};
Point(p22r+12) = {-23.2, -halfWidth, 35};
//Point(p22r+13) = {-0, -halfWidth, 35};

l33r = newl;
//Line (l33r+1) = {p22r+1, p22r+2} ;
Line (l33r+2) = {p22r+2, p22r+3} ;
Line (l33r+3) = {p22r+3, p22r+4} ;
Line (l33r+4) = {p22r+4, p22r+5} ;
Line (l33r+5) = {p22r+5, p22r+6} ;
Line (l33r+6) = {p22r+6, p22r+7} ;
Line (l33r+7) = {p22r+7, p22r+8} ;
Line (l33r+8) = {p22r+8, p22r+9} ;
Line (l33r+9) = {p22r+9, p22r+10} ;
Line (l33r+10) = {p22r+10, p22r+11} ;
Line (l33r+11) = {p22r+11, p22r+12} ;
//Line (l33+12) = {p22r+12, p22r+13} ;
//Line (l33+13) = {p22r+13, p22r+1} ;

Line (l33r+14) = {p22r+2, p2r+1} ;
Line (l33r+15) = {p22r+12, p2r+13} ;


lm3 = newl;
Line (lm3+1) = {p2+1,p2r+1} ;
Line (lm3+2) = {p2+2,p2r+2} ;
Line (lm3+3) = {p2+3,p2r+3} ;
Line (lm3+4) = {p2+4,p2r+4} ;
Line (lm3+5) = {p2+5,p2r+5} ;
Line (lm3+6) = {p2+6,p2r+6} ;
Line (lm3+7) = {p2+7,p2r+7} ;
Line (lm3+8) = {p2+8,p2r+8} ;
Line (lm3+9) = {p2+9,p2r+9} ;
Line (lm3+10) = {p2+10,p2r+10} ;
Line (lm3+11) = {p2+11,p2r+11} ;
Line (lm3+12) = {p2+12,p2r+12} ;
Line (lm3+13) = {p2+13,p2r+13} ;

lm4 = newl;
//Line (lm4+1) = {p22+1,p22r+1} ;
Line (lm4+2) = {p22+2,p22r+2} ;
Line (lm4+3) = {p22+3,p22r+3} ;
Line (lm4+4) = {p22+4,p22r+4} ;
Line (lm4+5) = {p22+5,p22r+5} ;
Line (lm4+6) = {p22+6,p22r+6} ;
Line (lm4+7) = {p22+7,p22r+7} ;
Line (lm4+8) = {p22+8,p22r+8} ;
Line (lm4+9) = {p22+9,p22r+9} ;
Line (lm4+10) = {p22+10,p22r+10} ;
Line (lm4+11) = {p22+11,p22r+11} ;
Line (lm4+12) = {p22+12,p22r+12} ;
//Line (lm4+13) = {p22+13,p22r+13} ;

//l4 = newll;
//Curve Loop(l4) = {l3+1, l3+2, l3+3, l3+4,l3+5,l3+6,l3+7,l3+8,l3+9,l3+10,l3+11,l3+12,l3+13};

//s2 = news;
//Plane Surface(s2) = {l4};


//coil
pcoil = newp;
Point (pcoil+1) = {9, 11, -24.6} ;
Point (pcoil+2) = {-9, 11, -24.6} ;
Point (pcoil+3) = {-9, -11, -24.6} ;
Point (pcoil+4) = {9, -11, -24.6} ;

lcoil = newl;
Line (lcoil+1) = {pcoil+1, pcoil+2} ;
Line (lcoil+2) = {pcoil+2, pcoil+3} ;
Line (lcoil+3) = {pcoil+3, pcoil+4} ;
Line (lcoil+4) = {pcoil+4, pcoil+1} ;


Point (pcoil+5) = {15, 17, -24.6} ;
Point (pcoil+6) = {-15, 17, -24.6} ;
Point (pcoil+7) = {-15, -17, -24.6} ;
Point (pcoil+8) = {15, -17, -24.6} ;
r=5;
Point (pcoil+9) = {15-r, 17-r, -24.6} ;
Point (pcoil+10) = {-15+r, 17-r, -24.6} ;
Point (pcoil+11) = {-15+r, -17+r, -24.6} ;
Point (pcoil+12) = {15-r, -17+r, -24.6} ;

Point (pcoil+13) = {15, 17-r, -24.6} ;
Point (pcoil+14) = {-15, 17-r, -24.6} ;
Point (pcoil+15) = {-15, -17+r, -24.6} ;
Point (pcoil+16) = {15, -17+r, -24.6} ;

Point (pcoil+17) = {15-r, 17, -24.6} ;
Point (pcoil+18) = {-15+r, 17, -24.6} ;
Point (pcoil+19) = {-15+r, -17, -24.6} ;
Point (pcoil+20) = {15-r, -17, -24.6} ;

//+
c1 = newc;
Circle(c1+1) = {pcoil+13, pcoil+9, pcoil+17};
//+
Circle(c1+2) = {pcoil+18, pcoil+10, pcoil+14};
//+
Circle(c1+3) = {pcoil+15, pcoil+11, pcoil+19};
//+
Circle(c1+4) = {pcoil+20, pcoil+12, pcoil+16};

lcoil = newl;
//+
Line(lcoil+5) = {pcoil+16, pcoil+13};
//+
Line(lcoil+6) = {pcoil+17, pcoil+18};
//+
Line(lcoil+7) = {pcoil+14, pcoil+15};
//+
Line(lcoil+8) = {pcoil+19, pcoil+20};

//coil line 2
coilheight = 49.2;
pcoil2 = newp;
Point (pcoil2+1) = {9, 11, -24.6+coilheight} ;
Point (pcoil2+2) = {-9, 11, -24.6+coilheight} ;
Point (pcoil2+3) = {-9, -11, -24.6+coilheight} ;
Point (pcoil2+4) = {9, -11, -24.6+coilheight} ;

lcoil2 = newl;
Line (lcoil2+1) = {pcoil2+1, pcoil2+2} ;
Line (lcoil2+2) = {pcoil2+2, pcoil2+3} ;
Line (lcoil2+3) = {pcoil2+3, pcoil2+4} ;
Line (lcoil2+4) = {pcoil2+4, pcoil2+1} ;


Point (pcoil2+5) = {15, 17, -24.6+coilheight} ;
Point (pcoil2+6) = {-15, 17, -24.6+coilheight} ;
Point (pcoil2+7) = {-15, -17, -24.6+coilheight} ;
Point (pcoil2+8) = {15, -17, -24.6+coilheight} ;
r=5;
Point (pcoil2+9) = {15-r, 17-r, -24.6+coilheight} ;
Point (pcoil2+10) = {-15+r, 17-r, -24.6+coilheight} ;
Point (pcoil2+11) = {-15+r, -17+r, -24.6+coilheight} ;
Point (pcoil2+12) = {15-r, -17+r, -24.6+coilheight} ;

Point (pcoil2+13) = {15, 17-r, -24.6+coilheight} ;
Point (pcoil2+14) = {-15, 17-r, -24.6+coilheight} ;
Point (pcoil2+15) = {-15, -17+r, -24.6+coilheight} ;
Point (pcoil2+16) = {15, -17+r, -24.6+coilheight} ;

Point (pcoil2+17) = {15-r, 17, -24.6+coilheight} ;
Point (pcoil2+18) = {-15+r, 17, -24.6+coilheight} ;
Point (pcoil2+19) = {-15+r, -17, -24.6+coilheight} ;
Point (pcoil2+20) = {15-r, -17, -24.6+coilheight} ;

//+
c11 = newc;
Circle(c11+1) = {pcoil2+13, pcoil2+9, pcoil2+17};
//+
Circle(c11+2) = {pcoil2+18, pcoil2+10, pcoil2+14};
//+
Circle(c11+3) = {pcoil2+15, pcoil2+11, pcoil2+19};
//+
Circle(c11+4) = {pcoil2+20, pcoil2+12, pcoil2+16};

lcoil2 = newl;
//+
Line(lcoil2+5) = {pcoil2+16, pcoil2+13};
//+
Line(lcoil2+6) = {pcoil2+17, pcoil2+18};
//+
Line(lcoil2+7) = {pcoil2+14, pcoil2+15};
//+
Line(lcoil2+8) = {pcoil2+19, pcoil2+20};

llinajie = newl;
Line(llinajie+1) = {pcoil+13,pcoil2+13};
Line(llinajie+2) = {pcoil+14,pcoil2+14};
Line(llinajie+3) = {pcoil+15,pcoil2+15};
Line(llinajie+4) = {pcoil+16,pcoil2+16};
Line(llinajie+5) = {pcoil+17,pcoil2+17};
Line(llinajie+6) = {pcoil+18,pcoil2+18};
Line(llinajie+7) = {pcoil+19,pcoil2+19};
Line(llinajie+8) = {pcoil+20,pcoil2+20};