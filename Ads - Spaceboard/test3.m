STARTMO = 1:6;
OMEGA = 0:30:150;

deg = pi/180;
R0 = 6378.14;

alt = 550; 
a = R0 + alt; e = 0; om = 0;
in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
u = 0*deg;

RES = zeros(numel(STARTMO),numel(OMEGA));

for i = 1:numel(STARTMO)
    gd0 = [2018 STARTMO(i) 1 12 0 0];
    fprintf('Start month %g -> ',STARTMO(i));
    for j = 1:numel(OMEGA)
        ini0 = [a e om OMEGA(j) in u];
        [TABLE, whichcity] = ad2('Silentmode',1,'starttime',gd0,'ini0',ini0);
        RES(i,j) = TABLE{end,end};
        fprintf('|');
    end
    fprintf('\n');
end

plot3(STARTMO, OMEGA, RES);