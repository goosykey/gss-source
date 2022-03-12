function sunlight( jd, showpoint )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    showpoint = true;
end

R0 = 6378.137;

rsun = astro.sun(jd);
rsun = rsun';

%rsun = rsun * 1.496e+8;

add3d.elem.plotCircle3D([0 0 0], rsun, R0, '--y');

light('position',rsun*1.496e+8, 'style','infinite');
lighting gouraud;

sundot = rsun*(R0+50)/sqrt(sum(rsun.^2,2));

if showpoint
    
    plot3(sundot(1),sundot(2),sundot(3),'py', 'markerfacecolor', 'y');
    
    text(sundot(1),sundot(2),sundot(3),'Sun','VerticalAlignment','bottom',...
        'HorizontalAlignment','center', 'color','y', 'FontSize',9);
end

end

