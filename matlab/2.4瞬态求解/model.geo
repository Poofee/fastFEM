/* Converted from AutoCad DXF file: model.DXF */
/* Tolerance 1e-20: 26 points / 29 curves */

u = 1; 
lc = 2e-3 ;
small = 1e-3;
medium = 10e-3;

Point (11) = {0.032 *u, 0.031 *u, 0 *u, lc} ;
Point (2) = {0.032 *u, 0.023 *u, 0 *u, lc} ;
Point (1) = {0.032 *u, -0.024 *u, 0 *u, lc} ;
Point (4) = {0.032 *u, -0.031 *u, 0 *u, lc} ;
Point (10) = {0.028 *u, 0.023 *u, 0 *u, lc} ;
Point (9) = {0.028 *u, -0.024 *u, 0 *u, lc} ;
Point (16) = {0.025 *u, 0.022 *u, 0 *u, lc} ;
Point (15) = {0.025 *u, -0.023 *u, 0 *u, lc} ;
Point (17) = {0.017 *u, 0.022 *u, 0 *u, lc} ;
Point (14) = {0.017 *u, -0.023 *u, 0 *u, lc} ;
Point (12) = {0.015 *u, 0.031 *u, 0 *u, small} ;
Point (13) = {0.015 *u, 0.023 *u, 0 *u, small} ;
Point (18) = {0.0145 *u, 0.0465 *u-disp, 0 *u, lc} ;
Point (19) = {0.0145 *u, -0.004 *u-disp, 0 *u, small} ;
Point (7) = {0.0145 *u, -0.011 *u, 0 *u, small} ;
Point (8) = {0.0145 *u, -0.024 *u, 0 *u, lc} ;
Point (24) = {0.008 *u, 0.0465 *u-disp, 0 *u, lc} ;
Point (23) = {0.008 *u, 0.0315 *u-disp, 0 *u, lc} ;
Point (20) = {0.005 *u, -0.0135 *u-disp, 0 *u, small} ;
Point (6) = {0.005 *u, -0.0205 *u, 0 *u, small} ;
Point (22) = {0.003 *u, 0.026 *u-disp, 0 *u, lc} ;
Point (21) = {0.003 *u, -0.0135 *u-disp, 0 *u, small} ;
Point (5) = {0.003 *u, -0.0205 *u, 0 *u, small} ;
Point (3) = {0.003 *u, -0.031 *u, 0 *u, lc} ;
Point (25) = {0 *u, 0.0465 *u-disp, 0 *u, lc} ;
Point (26) = {0 *u, -0.0135 *u-disp, 0 *u, lc} ;

Line (1) = {1, 2} ;
Line (2) = {3, 4} ;
Line (3) = {4, 1} ;
Line (4) = {5, 6} ;
Line (5) = {7, 8} ;
Line (6) = {8, 9} ;
Line (7) = {6, 7} ;
Line (8) = {10, 2} ;
Line (9) = {5, 3} ;
Line (10) = {10, 9} ;
Line (11) = {9, 1} ;
Line (12) = {2, 11} ;
Line (13) = {11, 12} ;
Line (14) = {12, 13} ;
Line (15) = {13, 10} ;
Line (16) = {14, 15} ;
Line (17) = {15, 16} ;
Line (18) = {16, 17} ;
Line (19) = {17, 14} ;
Line (20) = {18, 19} ;
Line (21) = {19, 20} ;
Line (22) = {20, 21} ;
Line (23) = {21, 22} ;
Line (24) = {22, 23} ;
Line (25) = {23, 24} ;
Line (26) = {24, 18} ;
Line (27) = {24, 25} ;
Line (28) = {25, 26} ;
Line (29) = {26, 21} ;

//+
Point(27) = {0, 0, 0, lc};
//+
Point(28) = {75e-3, 0-disp, 0, 10e-3};
//+
Point(29) = {0, 75e-3, 0, 10e-3};
//+
Point(30) = {0, -75e-3, 0, 10e-3};
//+
Point(31) = {0, -110e-3, 0, 10e-3};
//+
Point(32) = {0, 110e-3, 0, 10e-3};
//+
Point(33) = {110e-3, 0, 0, 10e-3};
//+
Circle(30) = {30, 27, 28};
//+
Circle(31) = {28, 27, 29};
//+
Circle(32) = {31, 27, 33};
//+
Circle(33) = {33, 27, 32};
//+
Line(34) = {31, 30};
//+
Line(35) = {30, 26};
//+
Line(36) = {25, 29};
//+
Line(37) = {29, 32};

//+
Curve Loop(1) = {2, 3, -11, -6, -5, -7, -4, 9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {16, 17, 18, 19};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 1, -8, 10};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {15, 8, 12, 13, 14};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {21, 22, 23, 24, 25, 26, 20};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {29, 23, 24, 25, 27, 28};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {30, 31, -36, -27, 26, 20, 21, 22, -29, -35};
//+
Curve Loop(8) = {2, 3, 1, 12, 13, 14, 15, 10, -6, -5, -7, -4, 9};
//+
Plane Surface(7) = {2, 7, 8};
//+
Curve Loop(9) = {32, 33, -37, -31, -30, -34};
//+
Plane Surface(8) = {9};
