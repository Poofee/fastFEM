%2020-04-05
%by Poofee
%初步的想法是，不断地改变反射过程中的导纳
% 如果把z1的改变放到最后，结果也很有意思
close all
V = 10;%volatge source
Z = 0.05;%the inner resist of source
Z0 = 1;%the guess value  1  0.01
Z1 = 0.005;
U = [0:0.001:6.82,6.823:0.000001:6.824,6.824:0.001:10];
I = arrayfun(@r,U);
% line([0,10],[0,0],'LineWidth',1,'Color','k')
hold on
% plot(U,I,'r-');%the nonlinear resist
hold on
% plot(U,(10 - U) / Z);%the source line
xlabel('U','FontSize',20,'FontName','Times New Roman');
ylabel('I','FontSize',20,'FontName','Times New Roman');
axis equal
axis off
%  axis([0 10 0 10]);
% plot([-1 11],[0 0],'k');
% plot([0 0],[-10 210],'k');
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
    ua = (V/Z -ib + ub/Z0)/(1/Z+1/Z0);
    Ur = ua - Ui;
    ia = (10-ua)/Z;
    %
%     [normx, normy] = coord2norm(gca, [ub,ua],[ib,ia]);
%     annotation('arrow', normx, normy);
    
    if textcount <10
       text(ua,ia,num2str(textcount)) ;
       textcount = textcount+1;
    end
    if textcount > 7
        plot(ua, ia,'or','MarkerSize',0.2,'markerfacecolor', 'black','MarkerEdgeColor', 'black');
        line([ub,ua],[ib,ia],'LineWidth',0.1,'Color','k','LineStyle',':')
    else
        plot(ua, ia,'or','MarkerSize',1.6,'markerfacecolor', 'black','MarkerEdgeColor', 'black');
        line([ub,ua],[ib,ia],'LineWidth',0.4,'Color','k','LineStyle',':')
    end
%     line([0,ub],[-2*Ui/Z0,ib],'LineWidth',1,'Color','r','LineStyle','--')

    drawnow
    pause(1)
%     saveas(gcf,[num2str(count) '.jpg']);
    count = count+1;
%     line([ua,ua],[ia,r(ua)],'linestyle',':')
%     plot(ua,r(ua),'or','MarkerSize',5,'markerfacecolor','blue')
%     pause(1)
%     saveas(gcf,[num2str(count) '.jpg']);
%     count = count+1;
%     arrow3([ub,ib],[ua,ia],'-k',1)
    %the reflect process
    y0 = r(ua)/ua;
    

    if count < 40
        Z1 = ua/ia;
    end
    ub = fzero(@(x)(r(x)-ia-1/Z1*(ua-x)),3)
    Ui = ub - Ur;
    ib = r(ub);
    
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
    
    if textcount <10
       text(ub,ib,num2str(textcount)) ;
       textcount = textcount+1;
    end
    if textcount > 7
        plot(ub,ib,'or','MarkerSize',0.2,'markerfacecolor','black','MarkerEdgeColor', 'black');
        line([ua,ub],[ia,ib],'LineWidth',0.1,'Color','k','LineStyle',':')
    else
        plot(ub,ib,'or','MarkerSize',1.6,'markerfacecolor','black','MarkerEdgeColor', 'black');
        line([ua,ub],[ia,ib],'LineWidth',0.4,'Color','k','LineStyle',':')
    end
%     line([0,ua],[2*Ur/Z0,ia],'LineWidth',1,'Color','r','LineStyle','--')

    drawnow
%     pause(0.1)
%     title(['Step: ',num2str(i)]);
%     saveas(gcf,[num2str(count) '.jpg']);
    
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

