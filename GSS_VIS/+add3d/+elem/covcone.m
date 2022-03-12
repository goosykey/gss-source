function covcone( sat, alpha, conecolor )
%PLOT COVERAGE CONES
%   all satellites have to be the same altitude for now

if nargin < 3
    conecolor = 'c';
end

R0 = 6378.137;
alt = sqrt(sum(sat.^2,2)) - R0;

boob = sind(alpha) * (R0 + alt)/R0;

if boob > 1
    boob = 1;
end

eps = acosd(boob);
beta = 90-eps-alpha;
Rcov = R0 * sind(beta);
bcov = R0 * cosd(beta);
acov = R0 + alt - bcov;

t = [0;Rcov/acov];
[X,Y,Z] = add3d.elem.cylinder_new(t);

X = X*acov + sat(1); Y = Y*acov + sat(2); Z = Z*acov + sat(3);

cone = surf(X*1e3,Y*1e3,Z*1e3,'FaceAlpha',0.1,'FaceColor',conecolor,'EdgeColor','none');

r0 = [0,0,1];
r1 = -sat/(R0+alt);

dir = cross(r0,r1);
ang = acosd(dot(r0,r1));

rotate(cone,dir,ang,sat*1e3);


end

