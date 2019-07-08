% 20190708 By Poofee
% ˲̬�������
% ģ��ΪFLUX����

% 20190707 by Poofee
% ����gmsh�Բ�ͬλ�Ƶķ���������ѹ�������Ƿ񲻱�

% ������
MAG_CIR = [1,3,4];% ��·
COIL = [2];% ��Ȧ
AIR = [7];% ����
MOBIL_AIR = [6];% �������˶��Ŀ���
CORE = [5];% ����
INFINITE = [8];% ����Զ����

mindisp = 0;% ��Сλ��
maxdisp = 6e-3;% ���λ��
step = 1e-3;

close all;

time = [0:1e-3:15e-3,...
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

% B-H����
H = 0:1:100000;
B = arrayfun(@getB,H);
figure
plot(H,B);
% ��·����
U = 28;
Rcoil = 3;
CoilTurn = 225;

% ��е����
m = 0.25;


cur_disp = 0;
% ����ʱ�䲽ѭ��
for t=1:length(time)
    disp('-----------------------------');
    disp([num2str(t),'-th step, time is ',num2str(time(t)),'s']);
    % ����
    disp('mesh of domain......');
    cmd = ['gmsh.exe -setnumber disp ', num2str(cur_disp),' -2 -format msh2 model.geo '];
    [status,cmdout] = system(cmd);    
    % pause(1);
    mesh = load_gmsh2('model.msh');
    disp(['mesh done. ',num2str(mesh.nbNod),' nodes, ',num2str(mesh.nbTriangles),' triangles elements.']);
    % ���˲̬�ų�
    disp('start solve transient magnetic field......');
    
    disp('start nonlinear Newton iteration......');
    
    disp('calculate electromagnetic force on moving body......');
    
    disp('calculate mechanical variables......');
    
    
end
