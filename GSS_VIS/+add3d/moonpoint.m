function moonpoint( jd )
%MOONPOINT Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.137;

rmoon = astro.moon(jd);

moondot = rmoon*(R0+50)/sqrt(sum(rmoon.^2,2));

triplet = [0.3364  0.5012  0.7973];


plot3(moondot(1),moondot(2),moondot(3),'^c', 'markerfacecolor', triplet);

text(moondot(1),moondot(2),moondot(3),'Moon','VerticalAlignment','bottom',...
    'HorizontalAlignment','center', 'color',triplet, 'FontSize',9);


end

