function [B]=getB(H)
mu0 = 4*pi*1e-7;
mur = 500;
Js = 1.9;


% if H < 100
%     B = 0.0628 / 100 * H;
% else
    B = mu0*H+2*Js/pi*atan(pi*(mur-1)*mu0*H/2/Js);
% end
if B > 3
    B = 3;
end
end