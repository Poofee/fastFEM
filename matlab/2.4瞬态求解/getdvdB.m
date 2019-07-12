function [dvdB] = getdvdB(B)
H = getH(B);
if(B < 1e-9)
   disp('B ÎªÁã');
   return
end
dvdB = 1/B/getdBdH(H)-H/B/B;
end