% 20190708 By Poofee
% ˲̬�������
% ģ��ΪFLUX���̣�
% ����ѹ��Ϊ0

% 20190707 by Poofee
% ����gmsh�Բ�ͬλ�Ƶķ���������ѹ�������Ƿ񲻱�

% �ο����ף�
% 1.K. Hollaus, J. Sch?berl, and H. Silm��Transient Finite Element Simulation of Non-Linear Eddy Current
% Problems with Biot-Savart-Field of Voltage-Driven Coils��

%% ������
MAG_CIR = [1,3,4];% ��·
COIL = 2;% ��Ȧ
AIR = 7;% ����
MOBILE_AIR = 6;% �������˶��Ŀ���
CORE = 5;% ����
INFINITE = 8;% ����Զ����

%% ��������
COMPRESSIBLE_PART = [AIR,INFINITE];% ��ѹ������
FIXED_PART = [MAG_CIR,COIL,MOBILE_AIR];% �����ƶ�����
MOBILE_PART = CORE;% ���ƶ�����
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
time = time';% ��ɢʱ��
displacement = zeros(length(time),1);% λ��
velocity = zeros(length(time),1);% �ٶ�
acceleration = zeros(length(time),1);% ���ٶ�
yForce = zeros(length(time),1);% �����
current = zeros(length(time),1);% ����
phicoil = zeros(length(time),1);% ��Ȧ�Ĵ�ͨ

%% B-H����
H = 0:1:100000;
B = arrayfun(@getB,H);
sigma = 1/(2e-7);% ohm*m

% figure(1)
% plot(H,B);
% xlabel('H(A/m)');
% ylabel('B(T)');
% title('���Ų��ϵ�B-H����');
% drawnow;
%% ��·����
U = 24;
Rcoil = 3;
CoilTurn = 225;
Scoil = 8e-3*45e-3;
tau = CoilTurn/Scoil;

%% ��е����
m = 0.25;% ��������
mindisp = 0;% ��Сλ��
maxdisp = 6e-3;% ���λ��

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
% ��������������
figure;
axError = gca;%ʹ�õ�ǰ�figures�ϵ�������
title(axError,'���');%�����������ϵı���
set(axError,'NextPlot','add');%�����ڸ�������ӻ�ͼ�����ǲ�ɾ�����л�ͼ
axError.YScale = 'log';%����Y��Ϊ��������
axError.XLim = [0 inf];%����x�᷶Χ����0��ʼ���Ҳ��Զ�����
axError.YLim = [0 inf];%����y�᷶Χ�ײ���0��ʼ���ϲ��Զ�����
axError.Position = [ 0.1   0.1   0.8    0.8];%������������figure�е�λ��
axError.DataAspectRatioMode = 'auto';%�����������Զ����¿̶ȱ���
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

disp('��ʼ���ƽ��......');
h_result=figure;
ax1 = subplot(2,3,1);set(ax1,'NextPlot','add');
plot1 = plot(ax1,time(1),displacement(1),'o-');
title(ax1,'λ��');
hold on
ax2 = subplot(2,3,2);set(ax2,'NextPlot','add');
plot2 = plot(ax2,time(1),velocity(1),'o-');
title(ax2,'�ٶ�');
ax3 = subplot(2,3,3);set(ax3,'NextPlot','add');
plot3 = plot(ax3,time(1),yForce(1),'o-');
title(ax3,'�������');
ax4 = subplot(2,3,4);set(ax4,'NextPlot','add');
plot4 = plot(ax4,time(1),current(1),'o-');
title(ax4,'����');
ax5 = subplot(2,3,5);set(ax5,'NextPlot','add');
plot5 = plot(ax5,time(1),phicoil(1),'o-');
title(ax5,'��ͨ');
ax6 = subplot(2,3,6);set(ax6,'NextPlot','add');
plot6 = plot(ax6,time(1),acceleration(1),'o-');
title(ax6,'���ٶ�');
resultAxes = [ax1,ax2,ax3,ax4,ax5,ax6];
fluxresult(ax1,ax2,ax3,ax4,ax5,ax6);
t=2;
set(plot4,'XData',time(1:t),'YData',current(1:t));
set(h_result,'Position',[0 0 1 0.85],'Units','normalized');
set(h_result,'Position',[0 0 1 0.85],'Units','normalized');
drawnow
%% ����ʱ�䲽ѭ��
for t=2:length(time)
    disp('-----------------------------');
    disp(['��ʼ�� ',num2str(t-1),' ������, ʱ��Ϊ ',num2str(time(t)),' ��']);
    % ����
    disp('��ʼ����......');
    cmd = ['gmsh.exe -setnumber disp ', num2str(-position_flux( t-1)),' -2 -format msh2 model.geo '];
    [status,cmdout] = system(cmd);    
    % pause(1);
    mesh_0 = mesh;
    mesh = load_gmsh2('model.msh');
    disp(['��������. �� ',num2str(mesh.nbNod),' ���ڵ�, ',num2str(mesh.nbTriangles),' �����ǵ�Ԫ.']);
    % ���˲̬�ų�
    disp('��ʼ���˲̬�ų�......');
    if t==2
        Ak1 = zeros(mesh.nbNod,1);
        NL = mesh.TRIANGLES(:,1:3);
        Domain = mesh.ELE_TAGS((mesh.nbElm-mesh.nbTriangles+1):end,2);
        FixNL = NL;
        mu = 4*pi*1e-7*ones(length(Domain),1);
        % ɾ���������֣�����ɾ��Ϊ�˲��ı�˳�򣬷������
        FixNLIndex = 1:1:length(Domain);
        for idomain=length(Domain):-1:1
            if isempty(find(FIXED_MESH==Domain(idomain), 1))
                FixNL(idomain,:) = [];
                FixNLIndex(idomain) = [];
            end
        end
        disp(['һ���� ',num2str(length(FixNL)),' ���̶�����'])
    end
    Ak = Ak1;
    FixNLk = FixNL;
    AREA_0 = AREA;
    mu_0 = mu;
    FixNLIndexk = FixNLIndex;
    [Ak1,FixNL,AREA,currentt,FixNLIndex,mu] = magsolve(t,mesh,time,Ak,FixNLk,current,AREA_0,FixNLIndexk,mu_0,axResult,plotErrorA);
    
    disp('��ʼ�����Ӧ����......');
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
    disp(['��Ȧ��ͨΪ ',num2str(coilphi),' ����Ϊ ',num2str(current(t)),' A']);
    
    disp('��ʼ���������ϵĵ������......');
    
    disp('��ʼ�����е����......');
    velocity(t) = (displacement(t-1) - displacement(t))/(time(t) - time(t-1));
    acceleration(t) = (velocity(t-1) - velocity(t))/(time(t) - time(t-1));
    
    % �������ݣ�����ͼ
    set(plot1,'XData',time(1:t),'YData',displacement(1:t));  
    set(plot2,'XData',time(1:t),'YData',velocity(1:t));  
    set(plot3,'XData',time(1:t),'YData',yForce(1:t));  
    set(plot4,'XData',time(1:t),'YData',current(1:t));  
    set(plot5,'XData',time(1:t),'YData',phicoil(1:t));  
    set(plot6,'XData',time(1:t),'YData',acceleration(1:t));  
    drawnow
end
save(datestr(now,'yyyymmddHHMM'))

