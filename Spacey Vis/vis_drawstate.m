function [ H ] = vis_drawstate( CART510, jd, pla,spp, cdata, markerstyle, markerfacecolor, orbitstyle, plotearth )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    pla = 17;
    spp = 30;
end

if nargin < 5
    cdata = imread('EM_hd.jpg');
end

if nargin < 6
    markerstyle = 'or';
    markerfacecolor = [0.75 0.75 0.75];
    orbitstyle = 'r';
    plotearth = true;
end

if plotearth
    [ H ] = vis_earthdraw(jd, 'axes', 'none', 'prime', 'none','quality','hd', 'cdata', cdata);
end
    
%R0 = 6378.137;

X = CART510(:,1)*1e03;
Y = CART510(:,2)*1e03;
Z = CART510(:,3)*1e03;

plot3(X,Y,Z,markerstyle, 'markerfacecolor', markerfacecolor); % bull


%% LET US DRAW THE FUCKING PLANES

t = 0:0.01:2*pi;
Rsp = norm(CART510(1,1:3));

for i = 1:pla
    sat = (i-1)*spp + 1;
    NORMAL = cross(CART510(sat,1:3), CART510(sat,4:6));
    NORMAL = NORMAL/norm(NORMAL); % normalizing the normal vector
    A = NORMAL(1); B = NORMAL(2); C = NORMAL(3);
    X = Rsp/sqrt(A^2+C^2) * (C*cos(t) - A*B*sin(t)/norm(NORMAL));
    Y = Rsp*sqrt(A^2+C^2)*sin(t)/norm(NORMAL);
    Z = -Rsp/sqrt(A^2+C^2) * (A*cos(t) + B*C*sin(t)/norm(NORMAL));
    plot3(X*1e3,Y*1e3,Z*1e3,orbitstyle, 'linewidth', 1);
end



end

