% 20190708 By Poofee
% 瞬态特性求解
% 模型为FLUX例程

% 20190707 by Poofee
% 测试gmsh对不同位移的分网，不可压缩区域是否不变

% 区域编号
MAG_CIR = [1,3,4];% 磁路
COIL = [2];% 线圈
AIR = [7];% 空气
MOBILE_AIR = [6];% 随衔铁运动的空气
CORE = [5];% 衔铁
INFINITE = [8];% 无穷远区域

% 物理区域
COMPRESSIBLE_PART = [AIR,INFINITE];% 可压缩区域
FIXED_PART = [MAG_CIR,COIL,MOBILE_AIR];% 不可移动区域
MOBILE_PART = [CORE];% 可移动区域
FIXED_MESH = [MAG_CIR,COIL,CORE];

close all;

time = [0,0:1e-3:15e-3,...
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
phicoil = zeros(length(time),1);% 线圈的磁通

% B-H曲线
H = 0:1:100000;
B = arrayfun(@getB,H);
sigma = 2e-7;% ohm*m

figure
plot(H,B);
xlabel('H(A/m)');
ylabel('B(T)');
title('铁磁材料的B-H曲线');
drawnow;
% 电路参数
U = 28;
Rcoil = 3;
CoilTurn = 225;
Scoil = 0.36e-3;
tau = CoilTurn/Scoil;

% 机械参数
m = 0.25;% 衔铁质量
mindisp = 0;% 最小位移
maxdisp = 6e-3;% 最大位移

cur_disp = 0;
Ak = 0;
Ak1 = 0;
FixNL = 0;
% 进行时间步循环
for t=2:length(time)
    disp('-----------------------------');
    disp(['开始第 ',num2str(t-1),' 步迭代, 时间为 ',num2str(time(t)),' 秒']);
    % 分网
    disp('开始分网......');
    cmd = ['gmsh.exe -setnumber disp ', num2str(cur_disp),' -2 -format msh2 model.geo '];
    [status,cmdout] = system(cmd);    
    % pause(1);
    mesh = load_gmsh2('model.msh');
    disp(['分网结束. 共 ',num2str(mesh.nbNod),' 个节点, ',num2str(mesh.nbTriangles),' 个三角单元.']);
    % 求解瞬态磁场
    disp('开始求解瞬态磁场......');
    Ak = Ak1;
    FixNLk = FixNL;
    [Ak1,FixNL] = magsolve(mesh,time(t)-time(t-1),Ak,t,FixNLk);
    
    disp('开始计算感应电流......');
    
    disp('开始计算衔铁上的电磁吸力......');
    
    disp('开始计算机械变量......');
    
    
end

disp('开始绘制结果......');
