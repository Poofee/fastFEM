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

time = [0,0,0.00100000000000000,0.00200000000000000,0.00300000000000000,0.00400000000000000,0.00500000000000000,0.00600000000000000,0.00700000000000000,0.00800000000000000,0.00900000000000000,0.0100000000000000,0.0110000000000000,0.0120000000000000,0.0130000000000000,0.0140000000000000,0.0150000000000000,0.0151000000000000,0.0152000000000000,0.0153000000000000,0.0154000000000000,0.0155000000000000,0.0165000000000000,0.0170000000000000,0.0171000000000000,0.0172000000000000,0.0173000000000000,0.0174000000000000,0.0175000000000000,0.0185000000000000,0.0195000000000000,0.0200000000000000,0.0250000000000000,0.0300000000000000,0.0350000000000000,0.0400000000000000,0.0450000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000];
% time = [0,0:1e-3:15e-3,...
%     (15e-3+1e-4):1e-4:15.5e-3,...
%     (15.5e-3+1e-3):1e-3:17e-3,...
%     (17e-3+1e-4):1e-4:17.5e-3,...
%     (17.5e-3+1e-3):1e-3:20e-3,...
%     (20e-3+5e-3):5e-3:50e-3,...
%     (50e-3+10e-3):10e-3:0.1];
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
sigma = 1/2e-7;% ohm*m

% figure
% plot(H,B);
% xlabel('H(A/m)');
% ylabel('B(T)');
% title('铁磁材料的B-H曲线');
% drawnow;
% 电路参数
U = 24;
Rcoil = 3;
CoilTurn = 225;
Scoil = 8e-3*45e-3;
tau = CoilTurn/Scoil;

% 机械参数
m = 0.25;% 衔铁质量
mindisp = 0;% 最小位移
maxdisp = 6e-3;% 最大位移

position_flux = [0,-5.24642361495153e-06,-2.68101213531427e-05,-7.76591015679292e-05,-0.000170930856806491,-0.000319240672326960,-0.000534628907106191,-0.000828518596095625,-0.00121214520184128,-0.00169682271387370,-0.00229491542355361,-0.00302092126264144,-0.00389287031998499,-0.00493527738685913,-0.00618519707416504,-0.00631293055523800,-0.00644349837490657,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000];
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
    cmd = ['gmsh.exe -setnumber disp ', num2str(-position_flux(t-1)),' -2 -format msh2 model.geo '];
    [status,cmdout] = system(cmd);    
    % pause(1);
    mesh = load_gmsh2('model.msh');
    disp(['分网结束. 共 ',num2str(mesh.nbNod),' 个节点, ',num2str(mesh.nbTriangles),' 个三角单元.']);
    % 求解瞬态磁场
    disp('开始求解瞬态磁场......');
    Ak = Ak1;
    FixNLk = FixNL;
    [Ak1,FixNL,AREA,currentt] = magsolve(t,mesh,time,Ak,t,FixNLk,current,phicoil(t-1));
    
    disp('开始计算感应电流......');
    coilphi = 0;
    for j=1:mesh.nbTriangles
        if(find(COIL == mesh.ELE_TAGS(mesh.nbElm-mesh.nbTriangles+j,2)))
            n = mesh.TRIANGLES(j,1:3);
            coilphi = coilphi + sum(Ak1(n))*AREA(j)/3;
        end
    end
    coilphi = coilphi * 2 * pi * tau;
    phicoil(t) = coilphi;
    current(t) = currentt;%1/3*(24 + (phicoil(t-1) - phicoil(t))/(time(t) - time(t-1)));
    disp(['线圈磁通为 ',num2str(coilphi),' 电流为 ',num2str(current(t)),' A']);
    
    disp('开始计算衔铁上的电磁吸力......');
    
    disp('开始计算机械变量......');
    
    
end
save(datestr(now,'yyyymmddHHMM'))

disp('开始绘制结果......');
h_result=figure;
subplot(2,3,1);
plot(time,displacement,'*-');
title('位移');
hold on
subplot(2,3,2);
plot(time,velocity,'*-');
title('速度');
subplot(2,3,3);
plot(time,yForce,'*-');
title('电磁吸力');
subplot(2,3,4);
plot(time,current,'*-');
title('电流');
subplot(2,3,5);
plot(time,phicoil,'*-');
title('磁通');
subplot(2,3,6);
plot(time,acceleration,'*-');
title('加速度');

fluxresult(h_result);

set(h_result,'Position',[0 0 1 0.85],'Units','normalized');
set(h_result,'Position',[0 0 1 0.85],'Units','normalized');
drawnow
