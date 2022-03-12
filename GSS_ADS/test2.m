clear;
clc;

deg = pi/180;

incl = 0:2.5:90;

MEANS0 = incl*0;
MEANS1 = incl*0;
MEANS2 = incl*0;
MEANS3 = incl*0;

for i = 1:numel(incl)
    fprintf('|');
    [~, ~, ~, ~, MEANS0(i)] = ad1('POI',[56 37],'silentmode',true, 'in', incl(i)*deg);
    [~, ~, ~, ~, MEANS1(i)] = ad1('POI',[41 -74],'silentmode',true, 'in', incl(i)*deg);
    [~, ~, ~, ~, MEANS2(i)] = ad1('POI',[34 -118],'silentmode',true, 'in', incl(i)*deg);
    [~, ~, ~, ~, MEANS3(i)] = ad1('POI',[22 114],'silentmode',true, 'in', incl(i)*deg);
end

fprintf('\n');

plot(incl,MEANS0);
grid;hold on;
plot(incl,MEANS1);
plot(incl,MEANS2);
plot(incl,MEANS3);
legend('Moscow | 56 deg','New York | 41 deg','Los Angeles | 34 deg','Hong Kong | 22 deg');