clear all;clc
 
mu = 1.26686535*10^17; % gravitational par. Jup
Rj = 69911000; % [m] Jupiter mean radii
Rp = 200000000; % [m] parigee
M0 = 1000;
Mdry = 500;
Isp = 3000;
Difference = 100;
v_inf = 21100 - 13100;
 
i = 0;
for dm = 1:0.00001:499
    i = i+1;
    rp = Rp - 200000*dm;
    eq_my = ( sqrt( (Isp*log(1000/(1000-dm)))^2 + 2*mu/rp)+ Isp*log((1000-dm)/500))^2 - 2*mu/rp;
    x(i) = dm;
    y(i) = eq_my;
    Difference = abs(eq_my - v_inf^2);
    if Difference <= 10
        dm
        break
    end
end
plot(x,y)
