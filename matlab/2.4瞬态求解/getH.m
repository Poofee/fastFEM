function [H] = getH(B)
error = 1e-6;

Hmin = 0;
Hmax = 1e5;

Hhalf = 0.5*(Hmin + Hmax);
while(abs(getB(Hhalf) - B) > error)
    if(getB(Hhalf) > B)
        Hmax = Hhalf;
    else
        Hmin = Hhalf;
    end
    Hhalf = 0.5*(Hmin + Hmax);
end
H = Hhalf;
end