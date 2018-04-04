i=0;
axis off
for z1=1:1:30
    for z2=1:1:30
        for z3=1:1:30
            step = findz(z1,z2,z3);
            if step<=2              
                [step,1/z1,1/z2,1/z3]
            end
            plot3(z1,z2,step,'.','MarkerSize',10);hold on
            pause(1e-3)
%             plot(i,step,'*');hold on
%             i=i+1;
        end        
    end
end
% for i=1:length(a)/100-2
%     plot(a(i*100:2:(i+2)*100));
%     pause(1)
% end