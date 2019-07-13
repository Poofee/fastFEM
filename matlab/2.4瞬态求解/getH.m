function [H] = getH(B)
if B > 3
    B = 3;
end
error = 1e-6;

Hmin = 0;
Hmax = 1e7;
% disp(num2str(B));
Hhalf = 0.5*(Hmin + Hmax);
while(abs(getB(Hhalf) - B) > error)
    if(getB(Hhalf) > B)
        Hmax = Hhalf;
    else
        Hmin = Hhalf;
    end
    if abs(Hmax-Hmin)<1e-3
        disp(['�Ҳ�����ӦB=',num2str(B),'��Hֵ']);
        return;
    end
    Hhalf = 0.5*(Hmin + Hmax);
end
H = Hhalf;
end