% 20190708 By Poofee
% 瞬态特性求解
% 模型为FLUX例程

% 20190707 by Poofee
% 测试gmsh对不同位移的分网，不可压缩区域是否不变

% 区域编号
MAG_CIR = [1,3,4];% 磁路
COIL = [2];% 线圈
AIR = [7];% 空气
MOBIL_AIR = [6];% 随衔铁运动的空气
CORE = [5];% 衔铁
INFINITE = [8];% 无穷远区域

mindisp = 0;% 最小位移
maxdisp = 6e-3;% 最大位移
step = 1e-3;

close all;

time = [0:1e-3:15e-3,...
    (15e-3+1e-4):1e-4:15.5e-3,...
    (15.5e-3+1e-3):1e-3:17e-3,...
    (17e-3+1e-4):1e-4:17.5e-3,...
    (17.5e-3+1e-3):1e-3:20e-3,...
    (20e-3+5e-3):5e-3:50e-3,...
    (50e-3+10e-3):10e-3:0.1];
time = time';% 离散时间
displacement = zeros(length(time),1);% 位移
velocity = zeros(length(time),1);% 速度
acceleration = zeros(length(time),1);% 加速度
yForce = zeros(length(time),1);% 电磁力
current = zeros(length(time),1);% 电流

% B-H曲线
H = 0:1:100000;
B = arrayfun(@getB,H);
figure
plot(H,B);
% 电路参数
U = 28;
Rcoil = 3;
CoilTurn = 225;

% 机械参数
m = 0.25;


cur_disp = 0;
% 进行时间步循环
for t=1:length(time)
    disp('-----------------------------');
    disp([num2str(t),'-th step, time is ',num2str(time(t)),'s']);
    % 分网
    disp('mesh of domain......');
    cmd = ['gmsh.exe -setnumber disp ', num2str(cur_disp),' -2 -format msh2 model.geo '];
    [status,cmdout] = system(cmd);    
    % pause(1);
    mesh = load_gmsh2('model.msh');
    disp(['mesh done. ',num2str(mesh.nbNod),' nodes, ',num2str(mesh.nbTriangles),' triangles elements.']);
    % 求解瞬态磁场
    disp('start solve transient magnetic field......');
    
    disp('start nonlinear Newton iteration......');
    
    disp('calculate electromagnetic force on moving body......');
    
    disp('calculate mechanical variables......');
    
    
end
