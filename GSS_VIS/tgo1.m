R0 = 6378.14;
deg = pi/180;

t = 0:15:86400*3;

R1 = R0 + 500;
in = astro.sunsyn(500)*deg;

Ra = R0+10000;
Rp = R0+200;

a1 = (Ra + Rp)/2;
e1 = (Ra-Rp)/(Ra+Rp);

ini0 = [R1, 0, 0, 0, in,0
R1, 0, 0, 10*deg, in, 10*deg
R1, 0, 0, 20*deg, in, 20*deg
R1, 0, 0, 30*deg, in, 30*deg];

ini1 = [a1, e1, 0, 0, 40*deg, 0];
ini2 = [a1, e1, 0, 90*deg, 50*deg, 0];

[eph, ephcart] = astro.J2circ_cube(ini0,t);

[ eph1 ] = astro.J2pert( ini1, t );
[R1,V1] = math.randv(eph1(:,1),eph1(:,2),eph1(:,5),eph1(:,4),eph1(:,3),eph1(:,6)-eph1(:,3));

[ eph2 ] = astro.J2pert( ini2, t );
[R2,V2] = math.randv(eph2(:,1),eph2(:,2),eph2(:,5),eph2(:,4),eph2(:,3),eph2(:,6)-eph2(:,3));

slice1 = permute([R1;V1],[3,2,1]);
slice2 = permute([R2;V2],[3,2,1]);

ephcart = [ephcart;slice1; slice2];

vis_movie(t, ephcart);