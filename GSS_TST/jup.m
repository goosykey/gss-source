clear; clc;

Mjup = 1.898e27;
G = 6.674e-20; %km

Isp = 300;
g0 = 10;
m0 = 1000;



m1 = 0:500;
m2 = 500-m1;

r = 200000 - 200*m1;

Vesc = sqrt(2*Mjup*G./r); %km/s

dvOberth = Isp * g0 * log((m0-m1)./(m0-500)) * 1e-3;

Vafter = dvOberth .* sqrt(1+2.*Vesc./dvOberth);

% plot(m1/500, Vafter);
% hold on; grid;
% plot(m1/500, Vesc);
% plot(m1/500, dvOberth);
% fprintf('Rmin = %g\n',r(end));
% legend('Vafter', 'Vesc', 'dVOberth');

[~, mini] = min(abs(Vafter-8.07));

mfin1 = m1(mini)
mfin2 = m2(mini)
perc1 = m1(mini)/(m1(mini)+m2(mini))
perc2 = m2(mini)/(m1(mini)+m2(mini))
ratio = m1(mini)/m2(mini)