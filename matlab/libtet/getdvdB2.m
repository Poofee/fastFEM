function [dvdB2,mu] = getdvdB2(Blist,Hlist,Bvalue)
% 20200420 by Poofee
% 根据B值，计算磁阻率对B2的偏导数，使用一阶插值
dvdB2 = 0;
mu = 1;
% 计算H
if Bvalue > Blist(end)
    Bvalue = Blist(end);
end
if Bvalue == 0
    return;
end
% 如果H太小的话，不能用绝对误差，因为H比绝对误差还小
error = 1e-6;
% 从BH表中查找H值
Hmin = 0;
Hmax = Hlist(end);
% disp(num2str(B));
Hhalf = 0.5*(Hmin + Hmax);
Btmp = interp1(Hlist,Blist,Hhalf);
while(abs(Btmp - Bvalue)/Bvalue > error)
    if(Btmp > Bvalue)
        Hmax = Hhalf;
    else
        Hmin = Hhalf;
    end
    Hhalf = 0.5*(Hmin + Hmax);
    Btmp = interp1(Hlist,Blist,Hhalf);
end
H = Hhalf;
mu = Bvalue/H;
% 计算nu对B2的偏导数，主要思路是将点附近线性化插值，用直线表示，然后代入求导。
for i = 1:length(Hlist)-1
    if Hlist(i) =< H && Hlist(i+1) >= H
        k = (Hlist(i)-Hlist(i+1))/((Blist(i)-Blist(i+1)));
        b = Hlist(i) - k*Blist(i);
        dvdB2 = -b*0.5/Bvalue^3;
    end
end
end