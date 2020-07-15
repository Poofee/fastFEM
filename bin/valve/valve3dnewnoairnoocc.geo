//SetFactory("OpenCASCADE");

maxsize = 20;
lc =2;

halfWidth = 8.05;
//down part
//left
p1 = newp;
Point(p1+1) = {0, -halfWidth, -34 , lc};
Point(p1+2) = {0, -halfWidth, -12.15, lc};
Point(p1+3) = {0.98, -halfWidth, -12.15, lc};
Point(p1+4) = {4.825, -halfWidth, 2.2, lc};
Point(p1+5) = {5.95, -halfWidth, 2.2, lc};
Point(p1+6) = {5.95, -halfWidth, -28, lc};
Point(p1+7) = {17.1, -halfWidth, -28, lc};
Point(p1+8) = {17.1, -halfWidth, 22.25, lc};
Point(p1+9) = {17.686, -halfWidth, 22.25, lc};
Point(p1+10) = {23, -halfWidth, 7.65, lc};
Point(p1+11) = {25.85, -halfWidth, 7.65, lc};
Point(p1+12) = {25.85, -halfWidth, 5.475, lc};
Point(p1+13) = {23, -halfWidth, 5.475, lc};
Point(p1+14) = {23, -halfWidth, 3.3, lc};
Point(p1+15) = {25.85, -halfWidth, 3.3, lc};
Point(p1+16) = {25.85, -halfWidth, -2.557, lc};
Point(p1+17) = {23.25, -halfWidth, -9.7, lc};
Point(p1+18) = {23.25, -halfWidth, -34, lc};

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
//Point(p11+1) = {0, -halfWidth, -34,lc};
//Point(p11+2) = {0, -halfWidth, -12.15,lc};
Point(p11+3) = {-0.98, -halfWidth, -12.15,lc};
Point(p11+4) = {-4.825, -halfWidth, 2.2,lc};
Point(p11+5) = {-5.95, -halfWidth, 2.2,lc};
Point(p11+6) = {-5.95, -halfWidth, -28,lc};
Point(p11+7) = {-17.1, -halfWidth, -28,lc};
Point(p11+8) = {-17.1, -halfWidth, 22.25,lc};
Point(p11+9) = {-17.686, -halfWidth, 22.25,lc};
Point(p11+10) = {-23, -halfWidth, 7.65,lc};
Point(p11+11) = {-25.85, -halfWidth, 7.65,lc};
Point(p11+12) = {-25.85, -halfWidth, 5.475,lc};
Point(p11+13) = {-23, -halfWidth, 5.475,lc};
Point(p11+14) = {-23, -halfWidth, 3.3,lc};
Point(p11+15) = {-25.85, -halfWidth, 3.3,lc};
Point(p11+16) = {-25.85, -halfWidth, -2.557,lc};
Point(p11+17) = {-23.25, -halfWidth, -9.7,lc};
Point(p11+18) = {-23.25, -halfWidth, -34,lc};

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
Point(p1r+1) = {0, halfWidth, -34,lc};
Point(p1r+2) = {0, halfWidth, -12.15,lc};
Point(p1r+3) = {0.98, halfWidth, -12.15,lc};
Point(p1r+4) = {4.825, halfWidth, 2.2,lc};
Point(p1r+5) = {5.95, halfWidth, 2.2,lc};
Point(p1r+6) = {5.95, halfWidth, -28,lc};
Point(p1r+7) = {17.1, halfWidth, -28,lc};
Point(p1r+8) = {17.1, halfWidth, 22.25,lc};
Point(p1r+9) = {17.686, halfWidth, 22.25,lc};
Point(p1r+10) = {23, halfWidth, 7.65,lc};
Point(p1r+11) = {25.85, halfWidth, 7.65,lc};
Point(p1r+12) = {25.85, halfWidth, 5.475,lc};
Point(p1r+13) = {23, halfWidth, 5.475,lc};
Point(p1r+14) = {23, halfWidth, 3.3,lc};
Point(p1r+15) = {25.85, halfWidth, 3.3,lc};
Point(p1r+16) = {25.85, halfWidth, -2.557,lc};
Point(p1r+17) = {23.25, halfWidth, -9.7,lc};
Point(p1r+18) = {23.25, halfWidth, -34,lc};

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
//Point(p11r+1) = {0, halfWidth, -34,lc};
//Point(p11r+2) = {0, halfWidth, -12.15,lc};
Point(p11r+3) = {-0.98, halfWidth, -12.15,lc};
Point(p11r+4) = {-4.825, halfWidth, 2.2,lc};
Point(p11r+5) = {-5.95, halfWidth, 2.2,lc};
Point(p11r+6) = {-5.95, halfWidth, -28,lc};
Point(p11r+7) = {-17.1, halfWidth, -28,lc};
Point(p11r+8) = {-17.1, halfWidth, 22.25,lc};
Point(p11r+9) = {-17.686, halfWidth, 22.25,lc};
Point(p11r+10) = {-23, halfWidth, 7.65,lc};
Point(p11r+11) = {-25.85, halfWidth, 7.65,lc};
Point(p11r+12) = {-25.85, halfWidth, 5.475,lc};
Point(p11r+13) = {-23, halfWidth, 5.475,lc};
Point(p11r+14) = {-23, halfWidth, 3.3,lc};
Point(p11r+15) = {-25.85, halfWidth, 3.3,lc};
Point(p11r+16) = {-25.85, halfWidth, -2.557,lc};
Point(p11r+17) = {-23.25, halfWidth, -9.7,lc};
Point(p11r+18) = {-23.25, halfWidth, -34,lc};

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
Point(p2+1) = {0, halfWidth, -11.15,lc};
Point(p2+2) = {0.595, halfWidth, -11.15,lc};
Point(p2+3) = {4.44, halfWidth, 3.2,lc};
Point(p2+4) = {5.5, halfWidth, 3.2,lc};
Point(p2+5) = {5.5, halfWidth, 30,lc};
Point(p2+6) = {18, halfWidth, 30,lc};
Point(p2+7) = {18, halfWidth, 23.625,lc};
Point(p2+8) = {23.45, halfWidth, 8.65,lc};
Point(p2+9) = {25.35, halfWidth, 8.65,lc};
Point(p2+10) = {25.35, halfWidth, 25,lc};
Point(p2+11) = {23.2, halfWidth, 25,lc};
Point(p2+12) = {23.2, halfWidth, 35,lc};
Point(p2+13) = {0, halfWidth, 35,lc};

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
//Point(p22+1) = {0, halfWidth, -11.15,lc};
Point(p22+2) = {-0.595, halfWidth, -11.15,lc};
Point(p22+3) = {-4.44, halfWidth, 3.2,lc};
Point(p22+4) = {-5.5, halfWidth, 3.2,lc};
Point(p22+5) = {-5.5, halfWidth, 30,lc};
Point(p22+6) = {-18, halfWidth, 30,lc};
Point(p22+7) = {-18, halfWidth, 23.625,lc};
Point(p22+8) = {-23.45, halfWidth, 8.65,lc};
Point(p22+9) = {-25.35, halfWidth, 8.65,lc};
Point(p22+10) = {-25.35, halfWidth, 25,lc};
Point(p22+11) = {-23.2, halfWidth, 25,lc};
Point(p22+12) = {-23.2, halfWidth, 35,lc};
//Point(p22+13) = {-0, halfWidth, 35,lc};

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
Point(p2r+1) = {0, -halfWidth, -11.15,lc};
Point(p2r+2) = {0.595, -halfWidth, -11.15,lc};
Point(p2r+3) = {4.44, -halfWidth, 3.2,lc};
Point(p2r+4) = {5.5, -halfWidth, 3.2,lc};
Point(p2r+5) = {5.5, -halfWidth, 30,lc};
Point(p2r+6) = {18, -halfWidth, 30,lc};
Point(p2r+7) = {18, -halfWidth, 23.625,lc};
Point(p2r+8) = {23.45, -halfWidth, 8.65,lc};
Point(p2r+9) = {25.35, -halfWidth, 8.65,lc};
Point(p2r+10) = {25.35, -halfWidth, 25,lc};
Point(p2r+11) = {23.2, -halfWidth, 25,lc};
Point(p2r+12) = {23.2, -halfWidth, 35,lc};
Point(p2r+13) = {0, -halfWidth, 35,lc};

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
//Point(p22r+1) = {0, -halfWidth, -11.15,lc};
Point(p22r+2) = {-0.595, -halfWidth, -11.15,lc};
Point(p22r+3) = {-4.44, -halfWidth, 3.2,lc};
Point(p22r+4) = {-5.5, -halfWidth, 3.2,lc};
Point(p22r+5) = {-5.5, -halfWidth, 30,lc};
Point(p22r+6) = {-18, -halfWidth, 30,lc};
Point(p22r+7) = {-18, -halfWidth, 23.625,lc};
Point(p22r+8) = {-23.45, -halfWidth, 8.65,lc};
Point(p22r+9) = {-25.35, -halfWidth, 8.65,lc};
Point(p22r+10) = {-25.35, -halfWidth, 25,lc};
Point(p22r+11) = {-23.2, -halfWidth, 25,lc};
Point(p22r+12) = {-23.2, -halfWidth, 35,lc};
//Point(p22r+13) = {-0, -halfWidth, 35,lc};

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
Point (pcoil+1) = {9, 11, -24.6,lc} ;
Point (pcoil+2) = {-9, 11, -24.6,lc} ;
Point (pcoil+3) = {-9, -11, -24.6,lc} ;
Point (pcoil+4) = {9, -11, -24.6,lc} ;

lcoil = newl;
Line (lcoil+1) = {pcoil+1, pcoil+2} ;
Line (lcoil+2) = {pcoil+2, pcoil+3} ;
Line (lcoil+3) = {pcoil+3, pcoil+4} ;
Line (lcoil+4) = {pcoil+4, pcoil+1} ;


Point (pcoil+5) = {15, 17, -24.6,lc} ;
Point (pcoil+6) = {-15, 17, -24.6,lc} ;
Point (pcoil+7) = {-15, -17, -24.6,lc} ;
Point (pcoil+8) = {15, -17, -24.6,lc} ;
r=5;
Point (pcoil+9) = {15-r, 17-r, -24.6,lc} ;
Point (pcoil+10) = {-15+r, 17-r, -24.6,lc} ;
Point (pcoil+11) = {-15+r, -17+r, -24.6,lc} ;
Point (pcoil+12) = {15-r, -17+r, -24.6,lc} ;

Point (pcoil+13) = {15, 17-r, -24.6,lc} ;
Point (pcoil+14) = {-15, 17-r, -24.6,lc} ;
Point (pcoil+15) = {-15, -17+r, -24.6,lc} ;
Point (pcoil+16) = {15, -17+r, -24.6,lc} ;

Point (pcoil+17) = {15-r, 17, -24.6,lc} ;
Point (pcoil+18) = {-15+r, 17, -24.6,lc} ;
Point (pcoil+19) = {-15+r, -17, -24.6,lc} ;
Point (pcoil+20) = {15-r, -17, -24.6,lc} ;

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
Point (pcoil2+1) = {9, 11, -24.6+coilheight,lc} ;
Point (pcoil2+2) = {-9, 11, -24.6+coilheight,lc} ;
Point (pcoil2+3) = {-9, -11, -24.6+coilheight,lc} ;
Point (pcoil2+4) = {9, -11, -24.6+coilheight,lc} ;

lcoil2 = newl;
Line (lcoil2+1) = {pcoil2+1, pcoil2+2} ;
Line (lcoil2+2) = {pcoil2+2, pcoil2+3} ;
Line (lcoil2+3) = {pcoil2+3, pcoil2+4} ;
Line (lcoil2+4) = {pcoil2+4, pcoil2+1} ;


Point (pcoil2+5) = {15, 17, -24.6+coilheight,lc} ;
Point (pcoil2+6) = {-15, 17, -24.6+coilheight,lc} ;
Point (pcoil2+7) = {-15, -17, -24.6+coilheight,lc} ;
Point (pcoil2+8) = {15, -17, -24.6+coilheight,lc} ;
r=5;
Point (pcoil2+9) = {15-r, 17-r, -24.6+coilheight,lc} ;
Point (pcoil2+10) = {-15+r, 17-r, -24.6+coilheight,lc} ;
Point (pcoil2+11) = {-15+r, -17+r, -24.6+coilheight,lc} ;
Point (pcoil2+12) = {15-r, -17+r, -24.6+coilheight,lc} ;

Point (pcoil2+13) = {15, 17-r, -24.6+coilheight,lc} ;
Point (pcoil2+14) = {-15, 17-r, -24.6+coilheight,lc} ;
Point (pcoil2+15) = {-15, -17+r, -24.6+coilheight,lc} ;
Point (pcoil2+16) = {15, -17+r, -24.6+coilheight,lc} ;

Point (pcoil2+17) = {15-r, 17, -24.6+coilheight,lc} ;
Point (pcoil2+18) = {-15+r, 17, -24.6+coilheight,lc} ;
Point (pcoil2+19) = {-15+r, -17, -24.6+coilheight,lc} ;
Point (pcoil2+20) = {15-r, -17, -24.6+coilheight,lc} ;

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

Line(llinajie+9) = {pcoil+1,pcoil2+1};
Line(llinajie+10) = {pcoil+2,pcoil2+2};
Line(llinajie+11) = {pcoil+3,pcoil2+3};
Line(llinajie+12) = {pcoil+4,pcoil2+4};

//+
Curve Loop(1) = {231, 241, 232, 238, 229, 239, 230, 240};
//+
Curve Loop(2) = {226, 227, 224, 225};
//+
Plane Surface(1) = {1, 2};


//+
Curve Loop(3) = {212, 222, 213, 219, 210, 220, 211, 221};
//+
Curve Loop(4) = {207, 208, 205, 206};
//+
Plane Surface(2) = {3, 4};
//+
Curve Loop(5) = {212, 249, -231, -245};
//+
Surface(3) = {5};
//+
Curve Loop(6) = {222, 250, -241, -249};
//+
Plane Surface(4) = {6};
//+
Curve Loop(7) = {213, 246, -232, -250};
//+
Surface(5) = {7};
//+
Curve Loop(8) = {219, 243, -238, -246};
//+
Plane Surface(6) = {8};
//+
Curve Loop(9) = {210, 247, -229, -243};
//+
Surface(7) = {9};
//+
Curve Loop(10) = {220, 248, -239, -247};
//+
Plane Surface(8) = {10};
//+
Curve Loop(11) = {211, 244, -230, -248};
//+
Surface(9) = {11};
//+
Curve Loop(12) = {221, 245, -240, -244};
//+
Plane Surface(10) = {12};
//+
Curve Loop(13) = {205, 252, -224, -251};
//+
Plane Surface(11) = {13};
//+
Curve Loop(14) = {206, 253, -225, -252};
//+
Plane Surface(12) = {14};
//+
Curve Loop(15) = {207, 254, -226, -253};
//+
Plane Surface(13) = {15};
//+
Curve Loop(16) = {208, 251, -227, -254};
//+
Plane Surface(14) = {16};



//up

//+
Curve Loop(17) = {169, 170, 171, 172, 176, -160, -159, -158, -157, -156, -155, -154, -153, -152, -151, -150, -149, -175, 163, 164, 165, 166, 167, 168};
//+
Plane Surface(15) = {17};
//+
Curve Loop(18) = {127, 128, 129, 130, 131, -147, -143, -142, -141, -140, -139, -138, -137, -136, -135, -134, 146, 120, 121, 122, 123, 124, 125, 126};
//+
Plane Surface(16) = {18};
//+
Curve Loop(19) = {170, -201, -141, 200};
//+
Plane Surface(17) = {19};
//+
Curve Loop(20) = {171, -202, -142, 201};
//+
Plane Surface(18) = {20};
//+
Curve Loop(21) = {172, -203, -143, 202};
//+
Plane Surface(19) = {21};
//+
Curve Loop(22) = {176, -160, -189, 131, -147, 203};
//+
Plane Surface(20) = {22};
//+
Curve Loop(23) = {159, -189, -130, 188};
//+
Plane Surface(21) = {23};
//+
Curve Loop(24) = {158, -188, -129, 187};
//+
Plane Surface(22) = {24};
//+
Curve Loop(25) = {157, -187, -128, 186};
//+
Plane Surface(23) = {25};
//+
Curve Loop(26) = {186, -156, -185, 127};
//+
Plane Surface(24) = {26};
//+
Curve Loop(27) = {155, -185, -126, 184};
//+
Plane Surface(25) = {27};
//+
Curve Loop(28) = {154, -184, -125, 183};
//+
Plane Surface(26) = {28};
//+
Curve Loop(29) = {153, -183, -124, 182};
//+
Plane Surface(27) = {29};
//+
Curve Loop(30) = {152, -182, -123, 181};
//+
Plane Surface(28) = {30};
//+
Curve Loop(31) = {151, -181, -122, 180};
//+
Plane Surface(29) = {31};
//+
Curve Loop(32) = {150, -180, -121, 179};
//+
Plane Surface(30) = {32};
//+
Curve Loop(33) = {179, -149, -175, -193, 146, 120};
//+
Plane Surface(31) = {33};
//+
Curve Loop(34) = {163, -194, -134, 193};
//+
Plane Surface(32) = {34};
//+
Curve Loop(35) = {164, -195, -135, 194};
//+
Plane Surface(33) = {35};
//+
Curve Loop(36) = {165, -196, -136, 195};
//+
Plane Surface(34) = {36};
//+
Curve Loop(37) = {166, -197, -137, 196};
//+
Plane Surface(35) = {37};
//+
Curve Loop(38) = {167, -198, -138, 197};
//+
Plane Surface(36) = {38};
//+
Curve Loop(39) = {168, -199, -139, 198};
//+
Plane Surface(37) = {39};
//+
Curve Loop(40) = {169, -200, -140, 199};
//+
Plane Surface(38) = {40};


//down

//+
Curve Loop(41) = {59, -79, -77, -76, -75, -74, -73, -72, -71, -70, -69, -68, -67, -66, -65, -64, -63, 80, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58};
//+
Plane Surface(39) = {41};
//+
Curve Loop(42) = {19, -39, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, 40, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
//+
Plane Surface(40) = {42};
//+
Curve Loop(43) = {19, -39, 118, 79, -59, -99};
//+
Plane Surface(41) = {43};
//+
Curve Loop(44) = {37, 118, -77, -117};
//+
Plane Surface(42) = {44};
//+
Curve Loop(45) = {36, 117, -76, -116};
//+
Plane Surface(43) = {45};
//+
Curve Loop(46) = {35, 116, -75, -115};
//+
Plane Surface(44) = {46};
//+
Curve Loop(47) = {34, 115, -74, -114};
//+
Plane Surface(45) = {47};
//+
Curve Loop(48) = {33, 114, -73, -113};
//+
Plane Surface(46) = {48};
//+
Curve Loop(49) = {32, 113, -72, -112};
//+
Plane Surface(47) = {49};
//+
Curve Loop(50) = {31, 112, -71, -111};
//+
Plane Surface(48) = {50};
//+
Curve Loop(51) = {30, 111, -70, -110};
//+
Plane Surface(49) = {51};
//+
Curve Loop(52) = {29, 110, -69, -109};
//+
Plane Surface(50) = {52};
//+
Curve Loop(53) = {109, -68, -108, 28};
//+
Plane Surface(51) = {53};
//+
Curve Loop(54) = {27, 108, -67, -107};
//+
Plane Surface(52) = {54};
//+
Curve Loop(55) = {26, 107, -66, -106};
//+
Plane Surface(53) = {55};
//+
Curve Loop(56) = {25, 106, -65, -105};
//+
Plane Surface(54) = {56};
//+
Curve Loop(57) = {24, 105, -64, -104};
//+
Plane Surface(55) = {57};
//+
Curve Loop(58) = {23, 104, -63, -103};
//+
Plane Surface(56) = {58};
//+
Curve Loop(59) = {40, 3, 84, -43, -80, -103};
//+
Plane Surface(57) = {59};
//+
Curve Loop(60) = {4, 85, -44, -84};
//+
Plane Surface(58) = {60};
//+
Curve Loop(61) = {5, 86, -45, -85};
//+
Plane Surface(59) = {61};
//+
Curve Loop(62) = {6, 87, -46, -86};
//+
Plane Surface(60) = {62};
//+
Curve Loop(63) = {7, 88, -47, -87};
//+
Plane Surface(61) = {63};
//+
Curve Loop(64) = {8, 89, -48, -88};
//+
Plane Surface(62) = {64};
//+
Curve Loop(65) = {9, 90, -49, -89};
//+
Plane Surface(63) = {65};
//+
Curve Loop(66) = {10, 91, -50, -90};
//+
Plane Surface(64) = {66};
//+
Curve Loop(67) = {11, 92, -51, -91};
//+
Plane Surface(65) = {67};
//+
Curve Loop(68) = {12, 93, -52, -92};
//+
Plane Surface(66) = {68};
//+
Curve Loop(69) = {13, 94, -53, -93};
//+
Plane Surface(67) = {69};
//+
Curve Loop(70) = {14, 95, -54, -94};
//+
Plane Surface(68) = {70};
//+
Curve Loop(71) = {15, 96, -55, -95};
//+
Plane Surface(69) = {71};
//+
Curve Loop(72) = {16, 97, -56, -96};
//+
Plane Surface(70) = {72};
//+
Curve Loop(73) = {17, 98, -57, -97};
//+
Plane Surface(71) = {73};
//+
Curve Loop(74) = {18, 99, -58, -98};
//+
Plane Surface(72) = {74};

//volume


//+
Surface Loop(1) = {15, 38, 17, 18, 19, 20, 21, 16, 24, 23, 22, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37};
//+
Surface Loop(2) = {8, 2, 3, 4, 5, 6, 7, 1, 9, 10, 12, 13, 14, 11};
//+
Volume(1) = {2};
//+
Volume(2) = {1};
//+
Surface Loop(3) = {71, 40, 41, 42, 39, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 72};
//+
Volume(3) = {3};

inf = 100;
// create air box around magnets
BoundingBox; // recompute model bounding box
cx = (General.MinX + General.MaxX) / 2;
cy = (General.MinY + General.MaxY) / 2;
cz = (General.MinZ + General.MaxZ) / 2;
lx = 2*inf + General.MaxX - General.MinX;
ly = 2*inf + General.MaxY - General.MinZ;
lz = 2*inf + General.MaxZ - General.MinZ;
p1 = newp; Point (p1) = {cx-lx/2, cy-ly/2, cz-lz/2, maxsize};
p2 = newp; Point (p2) = {cx+lx/2, cy-ly/2, cz-lz/2, maxsize};
l1 = newl; Line(l1) = {p1, p2};
e1[] = Extrude {0, ly, 0} { Line{l1}; };
e2[] = Extrude {0, 0, lz} { Surface{e1[1]}; };
Delete { Volume{e2[1]}; }//delete the box before create new one
ss[] = {e1[1],e2[0],e2[2],e2[3],e2[4],e2[5]};
sl1 = newsl; Surface Loop(sl1) = {ss[]};
vv[] = {sl1};

//skin[] = CombinedBoundary{ Volume{1}; };
//sl = newsl; Surface Loop(sl) = {skin[]};
//vv[] += sl;
//skin[] = CombinedBoundary{ Volume{2}; };
//sl = newsl; Surface Loop(sl) = skin[];
//vv[] += sl;
//skin[] = CombinedBoundary{ Volume{3}; };
//sl = newsl; Surface Loop(sl) = skin[];
//vv[] += sl;

v1 = newv; 
//Volume(v1) = {vv[]};

//+
Volume(284) = { 282,1, 2, 3};


