function covcircle( sat, alpha, formatstring )
%COVCIRCLE Summary of this function goes here
%   Detailed explanation goes here

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

magnitude = sqrt(sum(sat.^2,2));

subsat = sat*R0./magnitude;
center = sat*bcov/magnitude;

normal = sat - subsat;

add3d.elem.plotCircle3D(center*1e3, normal, Rcov*1e3, formatstring)

end

