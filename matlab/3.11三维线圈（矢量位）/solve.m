% 20191025 by Poofee
% ����ʸ���ⵥԪ������ά���ų��������
% ���ģ��Ϊ��Ȧ�����Ĵų�

% ��Ȧ��һЩ������
% ��Ȧ��ѹ��25 V
% ��Ȧ������1000
% ��Ȧ�ھ���16 mm
% ��Ȧ�⾶��24 mm
% ��Ȧ�߶ȣ�46 mm
close all
tic
%% ����
fprintf('��ʼ����......\n');
cmd = ['gmsh.exe  -3 -format msh2 coil.geo'];
[status,cmdout] = system(cmd);
mesh = load_gmsh2('coil.msh');
fprintf('��������. �� %d ����Ԫ. ���� %d ���ڵ�, %d ����, %d ��������, %d ��������\n',...
    mesh.nbElm,mesh.nbPoints,mesh.nbLines,mesh.nbTriangles,mesh.nbTets);
toc

%% ����
% DisplayNodes(mesh.POS);

% DisplayElements(mesh.TETS(:,1:4),mesh.POS,mesh.nbTets);

% DisplayEdgeSB(100,mesh.TETS(:,1:4),mesh.POS);

%% �������
Ucoil = 25;%��Ȧ��ѹ
Ncoil = 1000;%��Ȧ����
Rc = 10;%��Ȧ����
Scoil = 0.46 * (0.24 - 0.16);%��Ȧ���
Js = Ncoil*Ucoil/Rc/Scoil;%��Ȧ�����ܶ�
mu0 = 4*pi*1e-7;%�����ŵ���

%% ���㵥Ԫ����
volume = zeros(mesh.nbTets,1);
for i=1:mesh.nbTets
    X = mesh.POS(mesh.TETS(i,1:4),1);
    Y = mesh.POS(mesh.TETS(i,1:4),2);
    Z = mesh.POS(mesh.TETS(1,1:4),3);
    tmp = [ones(4,1),X,Y,Z];
    a = det(tmp)/6;
%     clf
%     if a < 0
% %         Lines = [1 2 3 4 4 1;...
% %              2 3 4 1 2 3];
% %         line(X(Lines),Y(Lines),Z(Lines),'Color',[0 0 0]);
%         fill3(X([1,2,3]),Y([1,2,3]),Z([1,2,3]),[0 0 1]);hold on
%         fill3(X([1,2,4]),Y([1,2,4]),Z([1,2,4]),[0 1 1]);
%         fill3(X([2,3,4]),Y([2,3,4]),Z([2,3,4]),[0 1 1]);
%         fill3(X([1,3,4]),Y([1,3,4]),Z([1,3,4]),[0 1 1]);
% %         text(X,Y,Z,{'1','2','3','4'},'Color',[1 0 0],'FontSize',24);
%         axis equal
%         drawnow
%     end
    volume(i) = a;
end



