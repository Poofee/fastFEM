%2016-11-05
%by Poofee
%
close all
V = 10;%volatge source
Z = 0.05;%the inner resist of source
Z0 = 0.008;%the guess value1
U = 0:0.0001:10;
I = arrayfun(@r,U);
% plot(U,I,'r-');%the nonlinear resist
hold on
% plot(U,(10 - U) / Z);%the source line
% xlabel('U','FontSize',20,'FontName','Times New Roman');
% ylabel('I','FontSize',20,'FontName','Times New Roman');
axis equal
axis off
%  axis([0 10 0 10]);
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
textcount=1;
for i = 1:10
    %the incident process
    %
    ua = (V/Z + 2*Ui/Z0)/(1/Z+1/Z0);
    Ur = ua - Ui;
    ia = ua/Z0-2*Ui/Z0;
    %
%     line([ub,ua],[ib,ia],'LineWidth',1,'Color','k')
%     plot(ua, ia,'or','MarkerSize',5,'markerfacecolor', 'blue')
    if textcount > 10
        plot(ua, ia,'ok','MarkerSize',0.8,'markerfacecolor', 'black','MarkerEdgeColor', 'black');
        line([ub,ua],[ib,ia],'LineWidth',0.1,'Color','k','LineStyle','-')
    else
        plot(ua, ia,'ok','MarkerSize',1.6,'markerfacecolor', 'black','MarkerEdgeColor', 'black');
        line([ub,ua],[ib,ia],'LineWidth',0.4,'Color','k','LineStyle','-')
    end
    if textcount <10
       text(ua,ia,num2str(textcount)) ;
       
    end
    textcount = textcount+1;
%     pause(1)
%     saveas(gcf,[num2str(count) '.jpg']);
    count = count+1;
    if textcount > 10
    line([ua,ua],[ia,r(ua)],'linestyle',':','LineWidth',0.1,'Color','k')
    plot(ua,r(ua),'ok','MarkerSize',0.8,'markerfacecolor','blue')
    else
    line([ua,ua],[ia,r(ua)],'linestyle',':','LineWidth',0.4,'Color','k')
    plot(ua,r(ua),'ok','MarkerSize',1.6,'markerfacecolor','blue')
    end
    if textcount <10
       text(ua,r(ua),num2str(textcount)) ;
       
    end
    textcount = textcount+1;
%     pause(1)
%     saveas(gcf,[num2str(count) '.jpg']);
    count = count+1;
%     arrow3([ub,ib],[ua,ia],'-k',1)
    %the reflect process
    y0 = r(ua)/ua;
    ub = 2*Ur/Z0/(y0 + 1/Z0);
    Ui = ub - Ur;
    ib = ub*y0;
    
%     arrow3([ua,ia],[ub,ib],'-k',1)
    if ub > ua
        if textcount > 10
        line([0,ub],[0,ib],'Color','red','LineWidth',0.1)
        else
            line([0,ub],[0,ib],'Color','red','LineWidth',0.4)
        end
%         pause(2)
%         saveas(gcf,[num2str(count) '.jpg']);
    count = count+1;
    else
        if textcount > 10
        line([0,ua],[0,r(ua)],'Color','red','LineWidth',0.1)
        else
            line([0,ua],[0,r(ua)],'Color','red','LineWidth',0.4)
        end
%         pause(2)
%         saveas(gcf,[num2str(count) '.jpg']);
    count = count+1;
    end
%     line([ua,ub],[ia,ib],'LineWidth',1,'Color','k')
%     plot(ub,ib,'or','MarkerSize',5,'markerfacecolor','blue');
     if textcount > 10
        plot(ub, ib,'ok','MarkerSize',0.8,'markerfacecolor', 'black','MarkerEdgeColor', 'black');
        line([ua,ub],[ia,ib],'LineWidth',0.1,'Color','k','LineStyle','-')
    else
        plot(ub, ib,'ok','MarkerSize',1.6,'markerfacecolor', 'black','MarkerEdgeColor', 'black');
        line([ua,ub],[ia,ib],'LineWidth',0.4,'Color','k','LineStyle','-')
    end
    if textcount <10
       text(ub,ib,num2str(textcount)) ;
       
    end
    textcount = textcount+1;
%     pause(1)
%     title(['Step: ',num2str(i)]);
%     saveas(gcf,[num2str(count) '.jpg']);
    count = count+1;
end

% ub/ib
%  for i=1:count-1
%     str = strcat(num2str(i), '.jpg');
%     A=imread(str);
%     [I,map]=rgb2ind(A,256);
%     if(i==1)
%         imwrite(I,map,'movefig.gif','DelayTime',0.1,'LoopCount',Inf)
%     else
%         imwrite(I,map,'movefig.gif','WriteMode','append','DelayTime',1)    
%     end
% end

