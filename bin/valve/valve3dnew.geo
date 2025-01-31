SetFactory("OpenCASCADE");

Point(1) = {0, 0, -34};
Point(2) = {0, 0, -12.15};
Point(3) = {0.98, 0, -12.15};
Point(4) = {4.825, 0, 2.2};
Point(5) = {5.95, 0, 2.2};
Point(6) = {5.95, 0, -28};
Point(7) = {17.1, 0, -28};
Point(8) = {17.1, 0, 22.25};
Point(9) = {17.686, 0, 22.25};
Point(10) = {23, 0, 7.65};
Point(11) = {25.85, 0, 7.65};
Point(12) = {25.85, 0, 5.475};
Point(13) = {23, 0, 5.475};
Point(14) = {23, 0, 3.3};
Point(15) = {25.85, 0, 3.3};
Point(16) = {25.85, 0, -2.557};
Point(17) = {23.25, 0, -9.7};
Point(18) = {23.25, 0, -34};
Line (1) = {1, 2} ;
Line (2) = {2, 3} ;
Line (3) = {3, 4} ;
Line (4) = {4, 5} ;
Line (5) = {5, 6} ;
Line (6) = {6, 7} ;
Line (7) = {7, 8} ;
Line (8) = {8, 9} ;
Line (9) = {9, 10} ;
Line (10) = {10, 11} ;
Line (11) = {11, 12} ;
Line (12) = {12, 13} ;
Line (13) = {13, 14} ;
Line (14) = {14, 15} ;
Line (15) = {15, 16} ;
Line (16) = {16, 17} ;
Line (17) = {17, 18} ;
Line (18) = {18, 1} ;
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
Plane Surface(1) = {1};

Point(19) = {0, 0, -11.15};
Point(20) = {0.595, 0, -11.15};
Point(21) = {4.44, 0, 3.2};
Point(22) = {5.5, 0, 3.2};
Point(23) = {5.5, 0, 30};
Point(24) = {18, 0, 30};
Point(25) = {18, 0, 23.625};
Point(26) = {23.45, 0, 8.65};
Point(27) = {25.35, 0, 8.65};
Point(28) = {25.35, 0, 25};
Point(29) = {23.2, 0, 25};
Point(30) = {23.2, 0, 35};
Point(31) = {0, 0, 35};
Line (19) = {19, 20} ;
Line (20) = {20, 21} ;
Line (21) = {21, 22} ;
Line (22) = {22, 23} ;
Line (23) = {23, 24} ;
Line (24) = {24, 25} ;
Line (25) = {25, 26} ;
Line (26) = {26, 27} ;
Line (27) = {27, 28} ;
Line (28) = {28, 29} ;
Line (29) = {29, 30} ;
Line (30) = {30, 31} ;
Line (31) = {31, 19} ;
Curve Loop(2) = {19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
Plane Surface(2) = {2};

Point (32) = {9, 11, -24.6} ;
Point (33) = {-9, 11, -24.6} ;
Point (34) = {-9, -11, -24.6} ;
Point (35) = {9, -11, -24.6} ;
Line (32) = {32, 33} ;
Line (33) = {33, 34} ;
Line (34) = {34, 35} ;
Line (35) = {35, 32} ;
Curve Loop(3) = {32, 33, 34, 35};

Point (36) = {15, 17, -24.6} ;
Point (37) = {-15, 17, -24.6} ;
Point (38) = {-15, -17, -24.6} ;
Point (39) = {15, -17, -24.6} ;
r=5;
Point (40) = {15-r, 17-r, -24.6} ;
Point (41) = {-15+r, 17-r, -24.6} ;
Point (42) = {-15+r, -17+r, -24.6} ;
Point (43) = {15-r, -17+r, -24.6} ;

Point (44) = {15, 17-r, -24.6} ;
Point (45) = {-15, 17-r, -24.6} ;
Point (46) = {-15, -17+r, -24.6} ;
Point (47) = {15, -17+r, -24.6} ;

Point (48) = {15-r, 17, -24.6} ;
Point (49) = {-15+r, 17, -24.6} ;
Point (50) = {-15+r, -17, -24.6} ;
Point (51) = {15-r, -17, -24.6} ;

//+
Circle(40) = {44, 40, 48};
//+
Circle(41) = {49, 41, 45};
//+
Circle(42) = {46, 42, 50};
//+
Circle(43) = {51, 43, 47};
//+
Line(44) = {47, 44};
//+
Line(45) = {48, 49};
//+
Line(46) = {45, 46};
//+
Line(47) = {50, 51};
//+
Curve Loop(4) = {47, 43, 44, 40, 45, 41, 46, 42};
//+
Curve Loop(5) = {34, 35, 32, 33};
//+
Plane Surface(3) = {4, 5};

Extrude {0, 8.05, 0} {Surface{1, 2}; }
Extrude {0, -8.05, 0} {Surface{1, 2}; }
Symmetry {1,0,0,0} { Duplicata{Volume{1:4};} }
BooleanUnion { Volume{1}; Delete;} { Volume{2:8}; Delete; }
Extrude {0, 0, 49.2} {Surface{3}; }

Box(4) = {-40, -40, -56, 80, 80, 112};
BooleanDifference{Volume{4}; Delete;}{Volume{1:3};}
