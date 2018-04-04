% 2017-10-30
% by Poofee
% ������ʹ�������η�������������Ų����ľ��ų�
% ʵ�ַ����ο����������������̵�ų���ֵ������P142-P146
clear all
close all
fname = ['mesh.mphtxt'];

[X,Y,NL,Domain]=readcomsoltri(fname);

% % �������еĵ�Ԫ
%  for i=1:length(NL)
%      plot(X(NL(i,[1 2 3 1])),Y(NL(i,[1 2 3 1])));
%      hold on;
%  end
%   axis equal

num_nodes = length(X);
num_elements = length(NL);

CE = zeros(3,3);      % CE --- �����洢ÿ����Ԫ��ϵ������
S = zeros(num_nodes,num_nodes);%ȫ�־���
F1 = zeros(num_nodes,1);

%AREA = zeros(num_elements,1);%���������
P = zeros(num_elements,3);
Q = zeros(num_elements,3);
R = zeros(num_elements,3);

XL = X(NL);
YL = Y(NL);

Q(:,1) = YL(:,2) - YL(:,3);
Q(:,2) = YL(:,3) - YL(:,1);
Q(:,3) = YL(:,1) - YL(:,2);

R(:,1) = XL(:,3) - XL(:,2);
R(:,2) = XL(:,1) - XL(:,3);
R(:,3) = XL(:,2) - XL(:,1);

P(:,1) = XL(:,2).*YL(:,3) - XL(:,3).*YL(:,2);
P(:,2) = XL(:,3).*YL(:,1) - XL(:,1).*YL(:,3);
P(:,3) = XL(:,1).*YL(:,2) - XL(:,2).*YL(:,1);

AREA = 0.5 * (Q(:,1).*R(:,2) - Q(:,2).*R(:,1));%���������

%��ģ����ֻ����������1Ϊ��������2Ϊ��Ȧ��
%����һ���������⡣
%�������Ų����ļ���
%B=mu0*mur*Hc+Br,Br=1.2,������֤x�Ḻ�������ŵ�mur=1.109
mu0 = 4*pi*1e-7;%�����ŵ���
Br = 1.2;
mur = 1.109;
Hc = Br /  mu0 / mur;
H = zeros(num_elements,1);%�������Ų����ļ�������
pmdomain = find(Domain == 2);%Ѱ����������ĵ�Ԫ
H(pmdomain) = Hc;%������������Ľ���������������Ϊ0



for i=1:num_elements
    for row = 1:3
        for col = 1:3
            CE(row,col) = (R(i,row)*R(i,col)+Q(i,row)*Q(i,col))/4/AREA(i)/mu0;
            S(NL(i,row),NL(i,col)) = S(NL(i,row),NL(i,col)) + CE(row,col);
        end
        k = row + 1;
        if row == 3
            k = 1;
        end
        F1(NL(i,row)) = F1(NL(i,row)) + cos(pi)*H(i)/2*(XL(i,row) - XL(i,k));
        F1(NL(i,k)) = F1(NL(i,k)) + cos(pi)*H(i)/2*(XL(i,row) - XL(i,k));
    end
end
%���ұ߽��
A = zeros(num_nodes,1);
freenodes = find(abs(X)~=1 & abs(Y) ~=1);
A(freenodes) = S(freenodes,freenodes)\F1(freenodes);

F = scatteredInterpolant(X,Y,A);
% figure
% F = scatteredInterpolant(X,Y,A(1:num_nodes));
tx = -1:1e-2:1;
ty = -1:1e-2:1;
[qx,qy] = meshgrid(tx,ty);
qz = F(qx,qy);
% mesh(qx,qy,qz);
% surf(qx,qy,qz);
figure
subplot(1,2,2);hold on
title('MATLAB');
contourf(qx,qy,qz,20);colorbar
axis equal

% ����COMSOL�������ļ���������ͼ����
fp = fopen('comsoldata.txt','r');
% ��ȡǰ��û�õ�9��
for i = 1:9
    fgets(fp);
end
% ��ȡ��������
comsoldata = fscanf(fp,'%lf %lf %lf\n',[ 3 num_nodes]);
comsoldata = comsoldata';
fclose(fp);

F = scatteredInterpolant(comsoldata(:,1),comsoldata(:,2),comsoldata(:,3));
% figure
% F = scatteredInterpolant(X,Y,A(1:num_nodes));
tx = -1:1e-2:1;
ty = -1:1e-2:1;
[qx,qy] = meshgrid(tx,ty);
qz = F(qx,qy);
% mesh(qx,qy,qz);
% surf(qx,qy,qz);
subplot(1,2,1);
contourf(qx,qy,qz,20);colorbar
title('COMSOL');
hold on
axis equal