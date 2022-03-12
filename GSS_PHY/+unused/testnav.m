function [ deltas3 ] = testnav( a,e,om,Om,in, delmax, epsmax, maxind )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

deltas3 = zeros((2*maxind-1)^2,8,8);
scrsz = get(groot,'ScreenSize');

for u = 0:pi/4:7*pi/4

[delv,epsv, deltas, deltaa ] = navdelta( a,e,om,Om,in,u, delmax, epsmax, maxind );


deltas3(:,:,u/(pi/4)+1) = deltas;

figure('Name','DELTA contours','NumberTitle','off', 'OuterPosition',[scrsz(3)/8,scrsz(4)/8,scrsz(3)*0.75,scrsz(4)*0.75]);
titlest = ['\Deltaa, m; u = ', num2str(u*180/pi)];
[C,h]=contourf(delv*1e3,epsv*1e3, deltaa');
title(titlest);
xlabel('\Deltaxyz, m')
ylabel('\Deltav, m/s')
clabel(C,h)
fnamest = ['deltaa',num2str(u*180/pi)];
print(fnamest, '-dpng'); 

end




end

