% 20190708 By Poofee
% ˲̬�������
% ģ��ΪFLUX����

% 20190707 by Poofee
% ����gmsh�Բ�ͬλ�Ƶķ���������ѹ�������Ƿ񲻱�

% ������
MAG_CIR = [1,3,4];% ��·
COIL = [2];% ��Ȧ
AIR = [7];% ����
MOBILE_AIR = [6];% �������˶��Ŀ���
CORE = [5];% ����
INFINITE = [8];% ����Զ����

% ��������
COMPRESSIBLE_PART = [AIR,INFINITE];% ��ѹ������
FIXED_PART = [MAG_CIR,COIL,MOBILE_AIR];% �����ƶ�����
MOBILE_PART = [CORE];% ���ƶ�����
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
time = time';% ��ɢʱ��
displacement = zeros(length(time),1);% λ��
velocity = zeros(length(time),1);% �ٶ�
acceleration = zeros(length(time),1);% ���ٶ�
yForce = zeros(length(time),1);% �����
current = zeros(length(time),1);% ����
phicoil = zeros(length(time),1);% ��Ȧ�Ĵ�ͨ

% B-H����
H = 0:1:100000;
B = arrayfun(@getB,H);
sigma = 1/2e-7;% ohm*m

% figure
% plot(H,B);
% xlabel('H(A/m)');
% ylabel('B(T)');
% title('���Ų��ϵ�B-H����');
% drawnow;
% ��·����
U = 24;
Rcoil = 3;
CoilTurn = 225;
Scoil = 8e-3*45e-3;
tau = CoilTurn/Scoil;

% ��е����
m = 0.25;% ��������
mindisp = 0;% ��Сλ��
maxdisp = 6e-3;% ���λ��

position_flux = [0,-5.24642361495153e-06,-2.68101213531427e-05,-7.76591015679292e-05,-0.000170930856806491,-0.000319240672326960,-0.000534628907106191,-0.000828518596095625,-0.00121214520184128,-0.00169682271387370,-0.00229491542355361,-0.00302092126264144,-0.00389287031998499,-0.00493527738685913,-0.00618519707416504,-0.00631293055523800,-0.00644349837490657,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000,-0.00650000000000000];
cur_disp = 0;
Ak = 0;
Ak1 = 0;
FixNL = 0;
% ����ʱ�䲽ѭ��
for t=2:length(time)
    disp('-----------------------------');
    disp(['��ʼ�� ',num2str(t-1),' ������, ʱ��Ϊ ',num2str(time(t)),' ��']);
    % ����
    disp('��ʼ����......');
    cmd = ['gmsh.exe -setnumber disp ', num2str(-position_flux(t-1)),' -2 -format msh2 model.geo '];
    [status,cmdout] = system(cmd);    
    % pause(1);
    mesh = load_gmsh2('model.msh');
    disp(['��������. �� ',num2str(mesh.nbNod),' ���ڵ�, ',num2str(mesh.nbTriangles),' �����ǵ�Ԫ.']);
    % ���˲̬�ų�
    disp('��ʼ���˲̬�ų�......');
    Ak = Ak1;
    FixNLk = FixNL;
    [Ak1,FixNL,AREA,currentt] = magsolve(t,mesh,time,Ak,t,FixNLk,current,phicoil(t-1));
    
    disp('��ʼ�����Ӧ����......');
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
    disp(['��Ȧ��ͨΪ ',num2str(coilphi),' ����Ϊ ',num2str(current(t)),' A']);
    
    disp('��ʼ���������ϵĵ������......');
    
    disp('��ʼ�����е����......');
    
    
end
save(datestr(now,'yyyymmddHHMM'))

disp('��ʼ���ƽ��......');
h_result=figure;
subplot(2,3,1);
plot(time,displacement,'*-');
title('λ��');
hold on
subplot(2,3,2);
plot(time,velocity,'*-');
title('�ٶ�');
subplot(2,3,3);
plot(time,yForce,'*-');
title('�������');
subplot(2,3,4);
plot(time,current,'*-');
title('����');
subplot(2,3,5);
plot(time,phicoil,'*-');
title('��ͨ');
subplot(2,3,6);
plot(time,acceleration,'*-');
title('���ٶ�');

fluxresult(h_result);

set(h_result,'Position',[0 0 1 0.85],'Units','normalized');
set(h_result,'Position',[0 0 1 0.85],'Units','normalized');
drawnow
