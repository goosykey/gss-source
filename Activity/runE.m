HR = 0; % HOURS OF UTC
TC = 3600; % SECONDS

vtau = 100;
dtau = 40;
stau = 400;

jd = jday(2017,3,20,HR,0,0);
deg = pi/180;

a = 6978.14; e = 0; om = 0; in = 97.792*deg;
Om = 20*deg; u = 45*deg;

dOm = 11.125*deg; du = 8*deg; ddu = 12*deg;
netdata = {dOm, du, ddu, vtau, dtau, stau};


[ SESDATA, SESCOORD, SESBEG ] = gettrafficearth( jd, TC, map, mapdata, BIGASSMASK );

clear HR TC;