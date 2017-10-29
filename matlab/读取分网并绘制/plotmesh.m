% 2017-10-29
% by Poofee
% ��ȡģ�͵����ǵ�Ԫ�����ļ���
% �����Ʒ���
clear all
close all

fname=['..\..\model\model_normal.mphtxt'];

[X,Y,NL,Domain]=readcomsoltri(fname);
 
% �������еĵ�Ԫ
 for i=1:length(NL)
     plot(X(NL(i,[1 2 3 1])),Y(NL(i,[1 2 3 1])));
     hold on;
 end
  axis equal
%  ���Ƶ������ڵĵ�Ԫ
figure
minDomain = min(Domain);
maxDomain = max(Domain);
for i=minDomain:maxDomain
    theDomain = find(Domain == i);
    for j=1:length(theDomain)
        fill(X(NL(theDomain(j),:)),Y(NL(theDomain(j),:)),[i*1/maxDomain i*1/maxDomain i*1/maxDomain]);
        hold on
    end
end
 axis equal