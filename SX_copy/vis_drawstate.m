function [ H ] = vis_drawstate( CART510, jd, pla,spp )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    pla = 17;
    spp = 30;
end

[ H ] = vis_earthdraw( jd);

R0 = 6378.137;

X = CART510(:,1)*1e03;
Y = CART510(:,2)*1e03;
Z = CART510(:,3)*1e03;

plot3(X,Y,Z,'or', 'markerfacecolor', 'c');

% 241 is the shit (-> 211 240 242 270 271 300)

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
    plot3(X*1e3,Y*1e3,Z*1e3,'c', 'linewidth', 1);
end


%% plot ECI axes

plot3([0 8e6], [0 0], [0 0], 'r', 'linewidth',2);
plot3([0 0], [0 8e6], [0 0], 'r', 'linewidth',2);
plot3([0 0], [0 0], [0 8e6], 'r', 'linewidth',2);
text(8e6,0,0,'X_E_C_I', 'color','r');
text(0,8e6,0,'Y_E_C_I', 'color','r');
text(0,0,8e6,'Z_E_C_I', 'color','r');


%% PLOT EQUATOR & PRIME MERID

t = 0:0.002:2*pi;

thetag = gstime(jd);

xeq = R0*cos(t);
xpr = R0*cos(t);
xgr = R0/sin(thetag) * (sin(thetag)*cos(thetag)*sin(t));
yeq = R0*sin(t);
ypr = 0*t;
ygr = R0*sin(thetag)*sin(t);
zeq = 0*t;
zpr = R0*sin(t);
zgr = -R0/sin(thetag) * (sin(thetag)*cos(t));

plot3(xeq*1e3, yeq*1e3, zeq*1e3, 'w', 'linewidth',2);
plot3(xpr*1e3, ypr*1e3, zpr*1e3, 'w', 'linewidth',2);
plot3(xgr*1e3, ygr*1e3, zgr*1e3, 'g', 'linewidth',2);





end

