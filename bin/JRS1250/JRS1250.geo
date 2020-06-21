GAP = 10E-3;
XINGCHENG = 1.75;
CS_DONGTIEXIN = -30;
CS_JINGTIEXIN = CS_DONGTIEXIN+XINGCHENG+25;
CS_XIAGAI = CS_JINGTIEXIN+27.5-0.05;
CS_SHANGGAI = CS_XIAGAI-39.22;
CS_XIANQUAN = CS_XIAGAI-4.98-2.5;
CS_DAOCIHUAN = CS_XIAGAI-4.98-2.5-17-2.5;
CS_YONGCI = CS_DAOCIHUAN;
PM_FX_X = 27;
PM_FX_Y = -8;

lc = 1e-2;
disp = 0;

maxsize = 20;
midsize = 1;
smallsize = 0.5;

Point(1) = {6.79/2, CS_DONGTIEXIN+38.17-GAP-20.98-11.14-disp, 0, maxsize};
Point(2) = {16.11/2, CS_DONGTIEXIN+38.17-GAP-20.98-11.14-disp, 0, maxsize};
Point(3) = {16.11/2, CS_DONGTIEXIN+38.17-GAP-20.98-disp, 0, maxsize};
Point(4) = {33.02/2, CS_DONGTIEXIN+38.17-GAP-20.98-disp, 0, maxsize};
Point(5) = {33.02/2, CS_DONGTIEXIN+38.17-GAP-20.98+2.46-disp, 0, maxsize};
Point(6) = {41.9/2, CS_DONGTIEXIN+38.17-GAP-20.98+2.46-disp, 0, lc};
Point(7) = {73.94/2, CS_XIAGAI+0-4.98+1.92, 0, maxsize};
Point(8) = {16.08/2, CS_DONGTIEXIN+38.17-GAP-6.48-disp, 0, maxsize};
Point(9) = {16.07/2, CS_JINGTIEXIN+22.34-11.05+6.78, 0, maxsize};
Point(10) = {41.84/2, CS_JINGTIEXIN+22.47-11.05, 0, lc};
Point(11) = {41.84/2, CS_JINGTIEXIN+22.47, 0, maxsize};
Point(12) = {7.99/2, CS_JINGTIEXIN+22.34-11.05+6.78, 0, maxsize};
Point(13) = {7.99/2, CS_XIAGAI+0-4.98, 0, maxsize};
Point(14) = {60.03/2, CS_XIAGAI+0-4.98, 0, maxsize};
Point(15) = {60.03/2, CS_XIAGAI+0-4.98+1.92, 0, maxsize};
Point(16) = {87.97/2, CS_XIAGAI+0-4.98+1.92, 0, maxsize};
Point(17) = {87.97/2, CS_XIAGAI+0, 0, maxsize};
Point(18) = {7.99/2, CS_XIAGAI+0, 0, maxsize};
Point(19) = {(41.92)/2, CS_DAOCIHUAN+0, 0, lc};
Point(20) = {16.07/2, CS_JINGTIEXIN+22.47-11.05, 0, lc};
Point(21) = {6.79/2, CS_DONGTIEXIN+38.17-GAP-6.48-disp, 0, maxsize};
Point(22) = {45.86/2, CS_DAOCIHUAN+0, 0, maxsize};
Point(23) = {45.86/2, CS_DAOCIHUAN+0-9.93, 0, maxsize};
Point(24) = {(41.92)/2, CS_DAOCIHUAN+0-9.93, 0, lc};
Point(25) = {16.08/2, CS_DONGTIEXIN+38.17-GAP-disp, 0, lc};
Point(26) = {24+6.015, CS_YONGCI+0, 0, maxsize};
Point(27) = {24+6.015, CS_YONGCI+0-9.93, 0, maxsize};
Point(29) = {20.01/2, CS_SHANGGAI+0, 0, maxsize};
Point(30) = {60.03/2, CS_SHANGGAI+0, 0, maxsize};
Point(31) = {60.03/2, CS_SHANGGAI-1.91, 0, maxsize};
Point(33) = {73.94/2, CS_SHANGGAI-8.06, 0, maxsize};
Point(34) = {20.01/2, CS_SHANGGAI-8.06, 0, maxsize};
Point(35) = {73.94/2, CS_SHANGGAI-1.91, 0, maxsize};
Point(36) = {41.9/2+2.5, CS_XIANQUAN+0, 0, maxsize};
Point(37) = {41.9/2+2.5+3, CS_XIANQUAN+0, 0, maxsize};
Point(38) = {41.9/2+2.5+3, CS_XIANQUAN+0-17, 0, maxsize};
Point(39) = {41.9/2+2.5, CS_XIANQUAN+0-17, 0, maxsize};
Point(40) = {0, -100, 0, maxsize};
Point(41) = {100, 0, 0, maxsize};
Point(42) = {0, 100, 0, maxsize};
Point(43) = {0, -120, 0, maxsize};
Point(44) = {120, 0, 0, maxsize};
Point(45) = {0, 120, 0, maxsize};
Point(46) = {34.9/2, CS_DONGTIEXIN+38.17-GAP-disp, 0, 10*lc};
Point(47) = {34.9/2, CS_DONGTIEXIN+38.17-GAP-1.91-disp, 0, 20*lc};
Point(48) = {41.9/2, CS_DONGTIEXIN+38.17-GAP-1.91-disp, 0, lc};
Point(49) = {0, 0, 0, maxsize};

Line (1) = {1, 2} ;
Line (2) = {2, 3} ;
Line (3) = {13, 12} ;
Line (4) = {4, 5} ;
Line (5) = {5, 6} ;
Line (6) = {1, 21} ;
Line (7) = {46, 47} ;
Line (8) = {21, 8} ;
Line (9) = {47, 48} ;
Line (10) = {10, 11} ;
Line (11) = {13, 11} ;
Line (12) = {9, 12} ;
Line (13) = {7, 35} ;
Line (14) = {11, 14} ;
Line (15) = {14, 15} ;
Line (16) = {15, 7} ;
Line (17) = {16, 17} ;
Line (18) = {17, 18} ;
Line (19) = {18, 13} ;
Line (20) = {19, 22} ;
Line (21) = {22, 26} ;
Line (22) = {23, 27} ;
Line (23) = {22, 23} ;
Line (24) = {23, 24} ;
Line (25) = {24, 19} ;
Line (26) = {35, 33} ;
Line (27) = {14, 26} ;
Line (28) = {26, 27} ;
Line (29) = {27, 30} ;
Line (30) = {9, 20} ;
Line (31) = {48, 6} ;
Line (32) = {30, 31} ;
Line (33) = {31, 35} ;
Line (34) = {20, 10} ;
Line (35) = {8, 25} ;
Line (36) = {33, 34} ;
Line (37) = {34, 29} ;
Line (38) = {25, 46} ;
Line (39) = {7, 16} ;
Line (40) = {36, 37} ;
Line (41) = {37, 38} ;
Line (42) = {38, 39} ;
Line (43) = {39, 36} ;
Circle (44) = {40, 49, 41} ;
Line (45) = {43, 40} ;
Line (46) = {45, 42} ;
Circle (47) = {42, 49, 41} ;
Line (48) = {42, 40} ;
Circle (49) = {43, 49, 44} ;
Circle (50) = {45, 49, 44} ;
Line (52) = {3, 4} ;
Line (53) = {29, 30} ;


//+
Curve Loop(1) = {1, 2, 52, 4, 5, -31, -9, -7, -38, -35, -8, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {42, 43, 40, 41};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {34, 10, -11, 3, -12, 30};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {53, 32, 33, 26, 36, 37};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {13, -33, -32, -29, -28, -27, 15, 16};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {21, 28, -22, -23};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {24, 25, 20, 23};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {14, 15, 16, 39, 17, 18, 19, 11};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {49, -50, 46, 47, -44, -45};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {44, -47, 48};
//+
Curve Loop(11) = {36, 37, 53, -29, -22, 24, 25, 20, 21, -27, -14, -10, -34, -30, 12, -3, -19, -18, -17, -39, 13, 26};
//+
Curve Loop(12) = {43, 40, 41, 42};
//+
Curve Loop(13) = {8, 35, 38, 7, 9, 31, -5, -4, -52, -2, -1, 6};
//+
Plane Surface(10) = {10, 11, 12, 13};


//
Mesh.CharacteristicLengthExtendFromBoundary = 1;
Mesh.CharacteristicLengthFactor = 1;
Mesh.CharacteristicLengthMin = 0;
Mesh.CharacteristicLengthMax = 5;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthFromPoints = 1;

//+
Characteristic Length {36, 37, 38, 39, 29, 34, 33, 35, 31, 30, 27, 23, 22, 26, 18, 17, 16, 7, 15, 14, 11, 13, 12, 9, 8, 3, 2, 1, 21, 4, 5} = midsize;
//Characteristic Length {20, 10} = smallsize;
