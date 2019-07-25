function [H] = getH(B)
if B > 3
    B = 3;
end
% 如果H太小的话，不能用绝对误差，因为H比绝对误差还小
error = 1e-6;

Hmin = 0;
Hmax = 1e7;
% disp(num2str(B));
Hhalf = 0.5*(Hmin + Hmax);
while(abs(getB(Hhalf) - B)/B > error)
    if(getB(Hhalf) > B)
        Hmax = Hhalf;
    else
        Hmin = Hhalf;
    end
%     if abs(Hmax-Hmin)<1e-3
%         disp(['找不到对应B=',num2str(B),'的H值']);
%         return;
%     end
    Hhalf = 0.5*(Hmin + Hmax);
end
H = Hhalf;
end