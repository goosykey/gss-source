function [ a,e,om,Om,in,u ] = xyz2kepler( x,y,z,vx,vy,vz )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = 398600.4418; % Гравитационный параметр земли (= G*M), км3/с2 <- РАЗМЕРНАЯ ВЕЛИЧИНА, КИЛОМЕТРЫ ТУТ
rxyz = [x,y,z]'; % радиус вектор
vxyz = [vx,vy,vz]'; % вектор скорости КА

Ixyz = cross(rxyz,vxyz); % какой-то там вектор моментов
I = norm(Ixyz); % его модуль

r = norm([x,y,z]); % модуль радиус-вектора
v = norm([vx,vy,vz]); % модуль вектора скорости

Esp = v^2/2 - mu/r; % потенциальная энергия

a = -mu/(2*Esp); % большая полуось орбиты
e = sqrt(1-I^2/(a*mu)); % эксцентриситет орбиты
p = a*(1-e^2); % фокальный параметр орбиты

in = acos(Ixyz(3)/I); % наклонение

Om = atan2(Ixyz(1),-Ixyz(2)); % долгота восходящего узла

u = atan2(z/sin(in), x*cos(Om)+y*sin(Om)); % аргумент широты

if abs(p-r) > 1E-04 % ЕСЛИ орбита не круговая
    nu = atan2(sqrt(p/mu)*dot(vxyz,rxyz),p-r); % истинная аномалия движения КА
    om = u - nu; % аргумент перицентра орбиты
else
    nu = 0; % ИНАЧЕ обе эти величины не имеют физического смысла,
    om = 0; % но это один из выходных параметров функции
end

end

