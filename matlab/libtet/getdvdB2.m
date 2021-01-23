function [dvdB2,mu] = getdvdB2(Blist,Hlist,Bvalue)
% 20200420 by Poofee
% ����Bֵ����������ʶ�B2��ƫ������ʹ��һ�ײ�ֵ
dvdB2 = 0;
mu = 1;
% ����H
if Bvalue > Blist(end)
    Bvalue = Blist(end);
end
if Bvalue == 0
    return;
end
% ���H̫С�Ļ��������þ�������ΪH�Ⱦ�����С
error = 1e-6;
% ��BH���в���Hֵ
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
% ����nu��B2��ƫ��������Ҫ˼·�ǽ��㸽�����Ի���ֵ����ֱ�߱�ʾ��Ȼ������󵼡�
for i = 1:length(Hlist)-1
    if Hlist(i) =< H && Hlist(i+1) >= H
        k = (Hlist(i)-Hlist(i+1))/((Blist(i)-Blist(i+1)));
        b = Hlist(i) - k*Blist(i);
        dvdB2 = -b*0.5/Bvalue^3;
    end
end
end