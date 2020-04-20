% 20190726 By Poofee
% 三维静磁场求解
% 模型为FLUX例程

clear all;
close all;

% 求解valve.feem定义的三维静磁场问题
tet_NL_Mag_Static_Newton('valve.feem');

