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

time = [0,0:1e-3:15e-3,...
    (15e-3+1e-4):1e-4:15.5e-3,...
    (15.5e-3+1e-3):1e-3:17e-3,...
    (17e-3+1e-4):1e-4:17.5e-3,...
    (17.5e-3+1e-3):1e-3:20e-3,...
    (20e-3+5e-3):5e-3:50e-3,...
    (50e-3+10e-3):10e-3:0.1];
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
sigma = 2e-7;% ohm*m

figure
plot(H,B);
xlabel('H(A/m)');
ylabel('B(T)');
title('���Ų��ϵ�B-H����');
drawnow;
% ��·����
U = 28;
Rcoil = 3;
CoilTurn = 225;
Scoil = 0.36e-3;
tau = CoilTurn/Scoil;

% ��е����
m = 0.25;% ��������
mindisp = 0;% ��Сλ��
maxdisp = 6e-3;% ���λ��

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
    cmd = ['gmsh.exe -setnumber disp ', num2str(cur_disp),' -2 -format msh2 model.geo '];
    [status,cmdout] = system(cmd);    
    % pause(1);
    mesh = load_gmsh2('model.msh');
    disp(['��������. �� ',num2str(mesh.nbNod),' ���ڵ�, ',num2str(mesh.nbTriangles),' �����ǵ�Ԫ.']);
    % ���˲̬�ų�
    disp('��ʼ���˲̬�ų�......');
    Ak = Ak1;
    FixNLk = FixNL;
    [Ak1,FixNL,AREA] = magsolve(mesh,time(t)-time(t-1),Ak,t,FixNLk,current(t-1),phicoil(t-1));
    
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
    disp(['��Ȧ��ͨΪ ',num2str(coilphi)]);
    
    disp('��ʼ���������ϵĵ������......');
    
    disp('��ʼ�����е����......');
    
    
end

disp('��ʼ���ƽ��......');
h=figure;
subplot(2,3,1);
plot(time,displacement);
title('λ��');
hold on
subplot(2,3,2);
plot(time,velocity);
title('�ٶ�');
subplot(2,3,3);
plot(time,yForce);
title('�������');
subplot(2,3,4);
plot(time,current);
title('����');
subplot(2,3,5);
plot(time,phicoil);
title('��ͨ');
subplot(2,3,6);
plot(time,acceleration);
title('���ٶ�');
drawnow
set(h,'Position',[0 0 1 0.85],'Units','normalized');
set(h,'Position',[0 0 1 0.85],'Units','normalized');