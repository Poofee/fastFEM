function [dBdH] = getdBdH(H)
mu0 = 4*pi*1e-7;
mur = 500;
Js = 1.9;

% B = mu0*H+2*Js/pi*atan(pi*(mur-1)*mu0*H/2/Js);
dBdH = mu0 + mu0*(mur-1)/(1+(pi*(mur-1)*mu0*H/2/Js)^2);
end