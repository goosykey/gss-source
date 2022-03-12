function [ output_args ] = megaplotXYZ(t, mm)

x = mm(:,1);
y = mm(:,2);
z = mm(:,3);
vx = mm(:,4);
vy = mm(:,5);
vz = mm(:,6);

L = length(x);
a = zeros(L,1);
e = zeros(L,1);
om = zeros(L,1);
Om = zeros(L,1);
in = zeros(L,1);
u = zeros(L,1);

for i = 1:L
    [ a(i),e(i),om(i),Om(i),in(i),u(i) ] = ...
        xyz2kepler1( x(i),y(i),z(i),vx(i),vy(i),vz(i) );
end


scrsz = get(groot,'ScreenSize');

figure('Name','Report: Semi-major axis','NumberTitle','off', 'OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
plot(t, a, 'r', 'LineWidth', 2); grid; % Построение графика и сетки
legend('a(t)', 0); % Легенда на графике
xlabel('time, sec');
ylabel('semi-major axis, km');


figure('Name','Report: RAAN & inclination','NumberTitle','off', 'OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
plot(t, [Om,in] .*180 ./pi, 'LineWidth', 2); grid;
legend('O(t)', 'i(t)', 0);
xlabel('time, sec');
ylabel('angles, deg');

figure('Name','Report: Eccentricity','NumberTitle','off', 'OuterPosition',[1 30 scrsz(3)/2 scrsz(4)/2]);
plot(t, e, 'm', 'LineWidth', 2); grid;
legend('e(t)', 0);
xlabel('time, sec');
ylabel('eccentricity');

figure('Name','Report: Arg. of perigee','NumberTitle','off', 'OuterPosition',[scrsz(3)/2 30 scrsz(3)/2 scrsz(4)/2]);
plot(t, om, 'c', 'LineWidth', 2); grid;
legend('o(t)', 0);
xlabel('time, sec');
ylabel('arg. of perigee, rad');

figure('Name','Report: Arg. of latitude','NumberTitle','off');
plot(t, u, 'LineWidth', 2); grid;
legend('u(t)', 0);
xlabel('time, sec');
ylabel('arg. of latitude, rad');

end

