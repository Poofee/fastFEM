% 2017-10-29
% by Poofee
% 读取模型的三角单元分网文件，
% 并绘制分网
clear all
close all

fname=['..\..\model\model_normal.mphtxt'];

[X,Y,NL,Domain]=readcomsoltri(fname);
 
% 绘制所有的单元
 for i=1:length(NL)
     plot(X(NL(i,[1 2 3 1])),Y(NL(i,[1 2 3 1])));
     hold on;
 end
  axis equal
%  绘制单个域内的单元
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