% 20190708 By Poofee
% 瞬态特性求解
% 模型为FLUX例程，
% 弹簧压力为0

% 20190707 by Poofee
% 测试gmsh对不同位移的分网，不可压缩区域是否不变

% 参考文献：
% 1.K. Hollaus, J. Sch?berl, and H. Silm《Transient Finite Element Simulation of Non-Linear Eddy Current
% Problems with Biot-Savart-Field of Voltage-Driven Coils》

%% 区域编号
MAG_CIR = [1,3,4];% 磁路
COIL = 2;% 线圈
AIR = 7;% 空气
MOBILE_AIR = 6;% 随衔铁运动的空气
CORE = 5;% 衔铁
INFINITE = 8;% 无穷远区域

%% 物理区域
COMPRESSIBLE_PART = [AIR,INFINITE];% 可压缩区域
FIXED_PART = [MAG_CIR,COIL,MOBILE_AIR];% 不可移动区域
MOBILE_PART = CORE;% 可移动区域
FIXED_MESH = [MAG_CIR,COIL,CORE];

close all;

time = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.0151,0.0152,0.0153,0.0154,0.0155,0.0165,0.017,0.0171,0.0172,0.0173,0.0174,0.0175,0.0185,0.0195,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.07,0.08,0.09,0.1];
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

%% B-H曲线
H = 0:1:100000;
B = arrayfun(@getB,H);
sigma = 1/(2e-7);% ohm*m

% figure(1)
% plot(H,B);
% xlabel('H(A/m)');
% ylabel('B(T)');
% title('铁磁材料的B-H曲线');
% drawnow;
%% 电路参数
U = 24;
Rcoil = 3;
CoilTurn = 225;
Scoil = 8e-3*45e-3;
tau = CoilTurn/Scoil;

%% 机械参数
m = 0.25;% 衔铁质量
mindisp = 0;% 最小位移
maxdisp = 6e-3;% 最大位移

position_flux = [0,-5.24642361495153e-06,-2.68101213531427e-05,-7.76591015679292e-05,-0.000170930856806491,-0.000319240672326960,-0.000534628907106191,-0.000828518596095625,-0.00121214520184128,-0.00169682271387370,-0.00229491542355361,-0.00302092126264144,-0.00389287031998499,-0.00493527738685913,-0.00618519707416504,-0.00631293055523800,-0.00644349837490657,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065,-0.0065];
cur_disp = 0;
Ak = 0;
Ak1 = 0;
FixNL = 0;
mu = 0;
mu_0 = 0;

mesh = 0;
mesh_0 = 0;

AREA = 0;
AREA_0 = 0;

FixNLIndex = 0;
% 设置误差输出曲线
figure;
axError = gca;%使用当前活动figures上的坐标轴
title(axError,'误差');%设置坐标轴上的标题
set(axError,'NextPlot','add');%允许在该轴上添加绘图，但是不删除现有绘图
axError.YScale = 'log';%设置Y轴为对数坐标
axError.XLim = [0 inf];%设置x轴范围左侧从0开始，右侧自动更新
axError.YLim = [0 inf];%设置y轴范围底部从0开始，上部自动更新
axError.Position = [ 0.1   0.1   0.8    0.8];%设置坐标轴在figure中的位置
axError.DataAspectRatioMode = 'auto';%设置坐标轴自动更新刻度比例
axError.PlotBoxAspectRatioMode = 'auto';

h = figure;
sizeScreen = get(0,'ScreenSize');
width = sizeScreen(3);
height = sizeScreen(4);
set(h,'Position',[(width-0.8*height*0.7)/2 48 0.8*height*0.7 0.8*height]);
axResult = gca;
axis(axResult,'equal');
set(axResult,'NextPlot','add');

plotErrorA = plot(axError,1,1,'-*');

current(1) = 0;%U / Rcoil;

disp('开始绘制结果......');
h_result=figure;
ax1 = subplot(2,3,1);set(ax1,'NextPlot','add');
plot1 = plot(ax1,time(1),displacement(1),'o-');
title(ax1,'位移');
hold on
ax2 = subplot(2,3,2);set(ax2,'NextPlot','add');
plot2 = plot(ax2,time(1),velocity(1),'o-');
title(ax2,'速度');
ax3 = subplot(2,3,3);set(ax3,'NextPlot','add');
plot3 = plot(ax3,time(1),yForce(1),'o-');
title(ax3,'电磁吸力');
ax4 = subplot(2,3,4);set(ax4,'NextPlot','add');
plot4 = plot(ax4,time(1),current(1),'o-');
title(ax4,'电流');
ax5 = subplot(2,3,5);set(ax5,'NextPlot','add');
plot5 = plot(ax5,time(1),phicoil(1),'o-');
title(ax5,'磁通');
ax6 = subplot(2,3,6);set(ax6,'NextPlot','add');
plot6 = plot(ax6,time(1),acceleration(1),'o-');
title(ax6,'加速度');
resultAxes = [ax1,ax2,ax3,ax4,ax5,ax6];
fluxresult(ax1,ax2,ax3,ax4,ax5,ax6);
t=2;
set(plot4,'XData',time(1:t),'YData',current(1:t));
set(h_result,'Position',[0 0 1 0.85],'Units','normalized');
set(h_result,'Position',[0 0 1 0.85],'Units','normalized');
drawnow
%% 进行时间步循环
for t=2:length(time)
    disp('-----------------------------');
    disp(['开始第 ',num2str(t-1),' 步迭代, 时间为 ',num2str(time(t)),' 秒']);
    % 分网
    disp('开始分网......');
    cmd = ['gmsh.exe -setnumber disp ', num2str(-position_flux( t-1)),' -2 -format msh2 model.geo '];
    [status,cmdout] = system(cmd);    
    % pause(1);
    mesh_0 = mesh;
    mesh = load_gmsh2('model.msh');
    disp(['分网结束. 共 ',num2str(mesh.nbNod),' 个节点, ',num2str(mesh.nbTriangles),' 个三角单元.']);
    % 求解瞬态磁场
    disp('开始求解瞬态磁场......');
    if t==2
        Ak1 = zeros(mesh.nbNod,1);
        NL = mesh.TRIANGLES(:,1:3);
        Domain = mesh.ELE_TAGS((mesh.nbElm-mesh.nbTriangles+1):end,2);
        FixNL = NL;
        mu = 4*pi*1e-7*ones(length(Domain),1);
        % 删掉空气部分，倒着删是为了不改变顺序，方便迭代
        FixNLIndex = 1:1:length(Domain);
        for idomain=length(Domain):-1:1
            if isempty(find(FIXED_MESH==Domain(idomain), 1))
                FixNL(idomain,:) = [];
                FixNLIndex(idomain) = [];
            end
        end
        disp(['一共有 ',num2str(length(FixNL)),' 个固定网格'])
    end
    Ak = Ak1;
    FixNLk = FixNL;
    AREA_0 = AREA;
    mu_0 = mu;
    FixNLIndexk = FixNLIndex;
    [Ak1,FixNL,AREA,currentt,FixNLIndex,mu] = magsolve(t,mesh,time,Ak,FixNLk,current,AREA_0,FixNLIndexk,mu_0,axResult,plotErrorA);
    
    disp('开始计算感应电流......');
    coilphi = 0;
    for j=1:mesh.nbTriangles
        if((COIL == mesh.ELE_TAGS(mesh.nbElm-mesh.nbTriangles+j,2)))
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
    velocity(t) = (displacement(t-1) - displacement(t))/(time(t) - time(t-1));
    acceleration(t) = (velocity(t-1) - velocity(t))/(time(t) - time(t-1));
    
    % 更新数据，并绘图
    set(plot1,'XData',time(1:t),'YData',displacement(1:t));  
    set(plot2,'XData',time(1:t),'YData',velocity(1:t));  
    set(plot3,'XData',time(1:t),'YData',yForce(1:t));  
    set(plot4,'XData',time(1:t),'YData',current(1:t));  
    set(plot5,'XData',time(1:t),'YData',phicoil(1:t));  
    set(plot6,'XData',time(1:t),'YData',acceleration(1:t));  
    drawnow
end
save(datestr(now,'yyyymmddHHMM'))

