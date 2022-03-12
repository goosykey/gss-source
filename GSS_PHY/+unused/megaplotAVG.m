function [ output_args ] = megaplotAVG(t, y)


p = y(:,1);
l1 = y(:,2);
l2 = y(:,3);
e = sqrt(l1.^2 + l2.^2);
a = p./(1-e.^2);
om = acos(l1./e);
scrsz = get(groot,'ScreenSize');

figure('Name','Report: Semi-major axis','NumberTitle','off', 'OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
plot(t, a, 'r', 'LineWidth', 2); grid; hold; % Построение графика и сетки
pa = polyfit(t,a,1);
aavg = pa(1)*t + pa(2);
plot(t, aavg, 'c', 'LineWidth', 2);
legend('a(t)', 0); % Легенда на графике
xlabel('time, sec');
ylabel('semi-major axis, km');


figure('Name','Report: RAAN & inclination','NumberTitle','off', 'OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
plot(t, y(:,4:5) .*180 ./pi, 'LineWidth', 2); grid; 
legend('O(t)', 'i(t)', 0);
xlabel('time, sec');
ylabel('angles, deg');

figure('Name','Report: Eccentricity','NumberTitle','off', 'OuterPosition',[1 30 scrsz(3)/2 scrsz(4)/2]);
plot(t, e, 'm', 'LineWidth', 2); grid; hold;
pe = polyfit(t,e,1);
eavg = pe(1)*t + pe(2);
plot(t, eavg, 'g', 'LineWidth', 2);
legend('e(t)', 0);
xlabel('time, sec');
ylabel('eccentricity');

figure('Name','Report: Arg. of perigee','NumberTitle','off', 'OuterPosition',[scrsz(3)/2 30 scrsz(3)/2 scrsz(4)/2]);
plot(t, om, 'c', 'LineWidth', 2); grid;
legend('o(t)', 0);
xlabel('time, sec');
ylabel('arg. of perigee, rad');

figure('Name','Report: Arg. of latitude','NumberTitle','off');
plot(t, y(:,6), 'LineWidth', 2); grid;
legend('u(t)', 0);
xlabel('time, sec');
ylabel('arg. of latitude, rad');

figure('Name','Report: LAMBDAS','NumberTitle','off');
plot(t, y(:,2:3), 'LineWidth', 2); grid;
legend('l1(t)','l2(t)', 0);
xlabel('time, sec');
ylabel('lambdas');

clc;

end

