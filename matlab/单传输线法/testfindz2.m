i=0;
close all
for z1=1:.01:3
    for z2=1:.01:3
            step = findz2(z1,z2);
            if step<=2              
                [step,1/z1,1/z2]
            end
%             plot3(z1,z2,step,'.','MarkerSize',10);hold on
%             pause(1e-3)
plot(i,step,'.');hold on
pause(1e-3);
i = i + 1;
    end
end