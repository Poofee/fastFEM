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
disp = 0;//0-2.2

maxsize = 5;

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
Point(46) = {34.9/2, CS_DONGTIEXIN+38.17-GAP-disp, 0, lc};
Point(47) = {34.9/2, CS_DONGTIEXIN+38.17-GAP-1.91-disp, 0, lc};
Point(48) = {41.9/2, CS_DONGTIEXIN+38.17-GAP-1.91-disp, 0, maxsize};
Point(49) = {0, 0, 0, maxsize};
Point(50) = {41.9/2, CS_DONGTIEXIN+38.17-GAP-1.91-disp-6.8, 0, lc};
Point(51) = {34/2, CS_SHANGGAI+0, 0, maxsize};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 50};
//+
Line(7) = {50, 48};
//+
Line(8) = {48, 47};
//+
Line(9) = {47, 46};
//+
Line(10) = {46, 25};
//+
Line(11) = {25, 8};
//+
Line(12) = {8, 21};
//+
Line(13) = {21, 1};
//+
Line(14) = {39, 38};
//+
Line(15) = {38, 37};
//+
Line(16) = {37, 36};
//+
Line(17) = {36, 39};
//+
Line(18) = {24, 23};
//+
Line(19) = {23, 22};
//+
Line(20) = {22, 19};
//+
Line(21) = {19, 24};
//+
Line(22) = {23, 27};
//+
Line(23) = {27, 26};
//+
Line(24) = {26, 22};
//+
Line(25) = {20, 10};
//+
Line(26) = {10, 11};
//+
Line(27) = {11, 13};
//+
Line(28) = {13, 12};
//+
Line(29) = {12, 9};
//+
Line(30) = {9, 20};
//+
Line(31) = {11, 14};
//+
Line(32) = {14, 15};
//+
Line(33) = {15, 7};
//+
Line(34) = {7, 16};
//+
Line(35) = {16, 17};
//+
Line(36) = {17, 18};
//+
Line(37) = {18, 13};
//+
Line(38) = {34, 33};
//+
Line(39) = {33, 35};
//+
Line(40) = {35, 7};
//+
Line(41) = {14, 26};
//+
Line(42) = {27, 30};
//+
Line(43) = {30, 31};
//+
Line(44) = {31, 35};
//+
Line(45) = {29, 51};
//+
Line(46) = {29, 34};
Circle (47) = {42, 49, 41} ;
Line (48) = {42, 40} ;
Circle (49) = {43, 49, 44} ;
Circle (50) = {45, 49, 44} ;
Circle (51) = {40, 49, 41} ;
Line (52) = {43, 40} ;
Line (53) = {45, 42} ;
Line(54) = {51, 30};

//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {14, 15, 16, 17};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {25, 26, 27, 28, 29, 30};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {31, 32, 33, 34, 35, 36, 37, -27};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {41, -23, 42, 43, 44, 40, -33, -32};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {22, 23, 24, -19};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {18, 19, 20, 21};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {38, 39, -44, -43, -54, -45, 46};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {49, -50, 53, 47, -51, -52};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {51, -47, 48};
//+
Curve Loop(11) = {38, 39, 40, 34, 35, 36, 37, 28, 29, 30, 25, 26, 31, 41, 24, 20, 21, 18, 22, 42, -54, -45, 46};
//+
Plane Surface(10) = {1, 2, 10, 11};

// control mesh size
Field[1] = Distance;
//Field[1].NodesList = {5};
Field[1].NNodesByEdge = 80;
Field[1].EdgesList = {10,25,6,21};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc*2;
Field[2].LcMax = 10;
Field[2].DistMin = lc*2;
Field[2].DistMax = 50;

Field[3] = Distance;
//Field[1].NodesList = {5};
Field[3].NNodesByEdge = 50;
Field[3].EdgesList = {3,45};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.2;
Field[4].LcMax = 10;
Field[4].DistMin = 1;
Field[4].DistMax = 50;

//Field[5] = Distance;
////Field[1].NodesList = {5};
//Field[5].NNodesByEdge = 25;
//Field[5].EdgesList = {44,47,48,49,50};
//
//Field[6] = Threshold;
//Field[6].IField = 5;
//Field[6].LcMin = 2;
//Field[6].LcMax = 5;
//Field[6].DistMin = 10;
//Field[6].DistMax = 20;

// Finally, let's use the minimum of all the fields as the background mesh field
Field[7] = Min;
Field[7].FieldsList = {2,4};

Background Field = 7;

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

//Mesh.CharacteristicLengthExtendFromBoundary = 1;
//Mesh.CharacteristicLengthFactor = 1;
//Mesh.CharacteristicLengthMin = 0.01;
//Mesh.CharacteristicLengthMax = 50;
//Mesh.CharacteristicLengthFromCurvature = 0;
//Mesh.CharacteristicLengthFromPoints = 1;





