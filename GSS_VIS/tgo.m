
% fig = figure;
% pX = [-1 1 0; 
%     0 0 0]';
% pY = [-1/3 -1/3 2/3;
%     -1/3 -1/3 2/3]';
% pZ = [0 0 0;
%     0 0.5 0]';
% 
% globe = patch(pX,pY,pZ,'red');
% globe.FaceAlpha = 0.5;
% 
% axis([-2 2 -2 2 -2 2])
% view(3)

[ fig, globe, misc ] = vis_earthdraw(2451545,'quality','low');


while(1)
    rotate([globe; misc],[0 0 1], 0.5);
    pause (0.04);
end