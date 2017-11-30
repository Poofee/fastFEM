%2017-11-09
%by Poofee
% 这个程序主要是验证新的思路
%可以通过调节Z0和Ui来观察收敛效果
close all
clear all
V = 10;%volatge source
Z = 0.05;%the inner resist of source
Z0 = 0.008;%the guess value  0.004 0.9
U = 0:0.0001:10;
I = arrayfun(@r,U);
% line([0,10],[0,0],'LineWidth',1,'Color','k')
hold on
plot(U,I,'k-');%the nonlinear resist
hold on
plot(U,(10 - U) / Z,'k-');%the source line
xlabel('U','FontSize',20,'FontName','Times New Roman');
ylabel('I','FontSize',20,'FontName','Times New Roman');
axis equal
axis off
axis([-1 11 -10 210]);
plot([-1 11],[0 0],'k');
plot([0 0],[-10 210],'k');
set(gca,'ytick','','yTickLabel','');
set(gca,'xtick','','xTickLabel','');
daspect('auto');
%10 = Ui + I * (Z + Z0)
Ui = 0;
Ur = 0;
ua = 0;
ub = 0;
ia = 0;
ib = 0;
%U = - Z0 * I + b
%b = U + Z0 * I
t = 0.1;
count=1;
for i = 1:20
    %the incident process
    %    
    ua = (V/Z + 2*Ui/Z0)/(1/Z+1/Z0);
    Ur = ua - Ui;
    ia = ua/Z0-2*Ui/Z0;
    %
    line([ub,ua],[ib,ia],'LineWidth',1,'Color','k')
    line([ua,ua],[ia,r(ua)],'LineWidth',1,'Color',[0.5 0.5 0.5],'LineStyle','--')%竖直线
    plot(ua, r(ua),'or','MarkerSize',3,'markerfacecolor', 'black','MarkerEdgeColor','k')
%     line([0,ub],[-2*Ui/Z0,ib],'LineWidth',1,'Color','r','LineStyle','--')
    plot(ua, ia,'or','MarkerSize',3,'markerfacecolor', 'black','MarkerEdgeColor','k')
%     pause(0.1)
    saveas(gcf,[num2str(count) '.jpg']);
    count = count+1;
%     line([ua,ua],[ia,r(ua)],'linestyle',':')
%     plot(ua,r(ua),'or','MarkerSize',5,'markerfacecolor','blue')
%     pause(1)
%     saveas(gcf,[num2str(count) '.jpg']);
%     count = count+1;
%     arrow3([ub,ib],[ua,ia],'-k',1)
    %the reflect process
    y0 = r(ua)/ua;
    heh = ua:0.01:(ua+3);
%     plot(heh,ua.*ua.*ua/5-ua*3/5*ua*ua+heh.*3/5*ua*ua,'b')
    ub = fzero(@(x)(1/Z0*x+x*dr(ua)+(r(ua)-ua*dr(ua))-2*Ur/Z0),3);
    Ui = ub - Ur;
    ib = r(ua)-ua*dr(ua)+ub*dr(ua);
    
%     arrow3([ua,ia],[ub,ib],'-k',1)
%     if ub > ua
%         line([0,ub],[0,ib],'Color',[1 0 0])
%         pause(2)
%         saveas(gcf,[num2str(count) '.jpg']);
%     count = count+1;
%     else
%         line([0,ua],[0,r(ua)],'Color',[1 0 0])
%         pause(2)
%         saveas(gcf,[num2str(count) '.jpg']);
%     count = count+1;
%     end
    line([ua,ub],[ia,ib],'LineWidth',1,'Color','k')
    line([ua,ub],[r(ua),ib],'LineWidth',0.5,'Color','k','LineStyle','--')
%     line([0,ua],[2*Ur/Z0,ia],'LineWidth',1,'Color','r','LineStyle','--')
    plot(ub,ib,'or','MarkerSize',3,'markerfacecolor','black','MarkerEdgeColor','k');
%     pause(0.1)
    title(['Step: ',num2str(i)]);
    saveas(gcf,[num2str(count) '.jpg']);
    count = count+1;
end

ub/ib
 for i=1:count-1
    str = strcat(num2str(i), '.jpg');
    A=imread(str);
    [I,map]=rgb2ind(A,256);
    if(i==1)
        imwrite(I,map,'movefig.gif','DelayTime',0.1,'LoopCount',Inf)
    else
        imwrite(I,map,'movefig.gif','WriteMode','append','DelayTime',1)    
    end
end

