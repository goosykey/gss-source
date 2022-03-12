
%% INI

a = int8(9);

b = 0;

R0 = 6378.137;

mu = 398600.44;

h = 600;

%% IMPL

a = R0 + h;

T = 2*pi*sqrt(a^3/mu); % ckefcyfekcf

if T < 7200
    disp ('LOL');
elseif T == 10000
    disp ('hehe');
else
    disp ('not lol');
end

%% CYCLE

a = 1:10;

for i = a
    b = b + 1;
end

while b > 0
    b = b - 1;
end