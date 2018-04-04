function [dvdB] = getdvdB(B)
% 2017-10-30
% by Poofee
% 铁芯BH曲线表
% 数据来自comsol材料库
% 该函数返回磁阻率对B的偏微分再除以B，采用线性插值
table = [
    0 0
    663.146 1
    1067.5 1.1
    1705.23 1.2
    2463.11 1.3
    3841.67 1.4
    5425.74 1.5
    7957.75 1.6
    12298.3 1.7
    20462.8 1.8
    32169.6 1.9
    61213.4 2.0
    111408 2.1
    175070 2.2
    261469 2.3
    318310 2.4];

Hdata = table(:,1);
Bdata = table(:,2);

if B < 1e-9
    dvdB = 0;
    return;
end
for i=1:length(Hdata)-1
    if B >= Bdata(i) && B <= Bdata(i+1)
        slope = (Hdata(i+1) - Hdata(i)) / (Bdata(i+1) - Bdata(i));
        H = Hdata(i) + slope *(B - Bdata(i));
        b = Hdata(i) - slope * Bdata(i);
        dvdB = -b / B / B / B;
        return;
    end
end
i = length(Hdata) - 1;
slope = (Hdata(i+1) - Hdata(i)) / (Bdata(i+1) - Bdata(i));
H = Hdata(i) + slope *(B - Bdata(i));
b = Hdata(i) - slope * Bdata(i);
dvdB = -b / B / B / B;
end

