gd0 = [2017 03 20 12 0 0]; % SET UP STARTTIME (UTC)
tau1 = 81.92; % specify type1 session duration
tau2 = 1.92;  % specify type2 session duration
tau3 = 20.48; % specify type3 session duration
averate1 = 388090; %byte/day
averate2 = 53860;  %byte/day
averate3 = 289540; %byte/day
br_m = 500;
br_f = 500;
jd = jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6)); % find starttime Julian Date
deg = pi/180;
a = 7178.14; e = 0; om = 0; in = sunsyn(800)*deg; % specify initial orbit parameters
Om = 0*deg; u = 45*deg;
dOm = 25.4286*deg; du = 11.4286*deg; ddu = 17.1429*deg;  % specify constellation parameters
netdata = {dOm, du, ddu, tau1, tau2, tau3, averate1, averate2, averate3, br_m, br_f};
R0 = 6378.14;

cfg = [17,35];

elevcrit = 45; % deg
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + 800));
Rcov = sqrt(2*R0^2*(1-cosd(90-alphacrit-elevcrit)));
M = repmat([3.333 1.734 1.6]*1e6,[3,1]);
[ map, mapdata ] = readmap( 'glp15ag.asc' );    % read population density map
[ map1, map2, map3 ] = gensubs( map, mapdata, M, 25, 0.5);

[ STATE168, CART168, LATLON168 ] = seedconstel( [a e om Om in u], jd, 17, 35 );

alldata = zeros(168,1);

for i = 1:168
    fprintf(' > %g',i);
    [ data1, data2, data3 ] = gettrafficsx( LATLON168(i,1), LATLON168(i,2), Rcov, -1, map1, map2, map3, mapdata, netdata );
    alldata(i) = sum([ data1, data2, data3 ]);
end

alldata_reduc = alldata ./ meanvis168(LATLON168(:,1));

latlonGS = [64.83, -147.75, 2 % Fairbanks, Alaska
    21.321605, -157.861396, 1 % Honolulu, Hawaii
    33.933216, -118.281481, 1 % Los Angeles, California
    43.017832, -89.385132, 1  % Madison, Wisconsin
    25.812291, -80.355975, 3  % Miami, Florida
    18.399474, -66.049074, 1  % San Juan, Puerto Rico
    9.859500, 8.887399, 1     % Jos, Nigeria
    48.388457, -4.510587, 1   % Brest, France
    41.634157, -4.738647, 1   % Valladolid, Spain
    46.034898, 14.538621, 1   % Ljubljana, Slovenija
    45.498342, 13.747929, 1   % Pomjan, Slovenija
    41.124640, 14.774235, 1   % Benevento, Italy
    55.749287, 37.624671, 3   % Moscow, Russia
    53.210687, 50.231718, 3   % Samara, Russia
    58.005860, 56.259677, 1   % Perm, Russia
    56.017824, 92.961860, 1   % Krasnoyarsk, Russia
    52.301306, 104.285353, 3  % Irkutsk, Russia
    59.557912, 150.815451, 2  % Magadan, Russia
    30.780302, 120.776480, 1  % Jiaxing, China
    ];

satarray = CART168(:,1:3);

[ hubnums ] = findhubs( latlonGS, jd, satarray, 15 );

aol = STATE168(:,6);
prob = 0.1;
[ SATCOOR, CONNout ] = generatenet( cfg, prob, aol);
[ DEPOTS, USAGE, isol, chhubs ] = ISCLoad( CONNout,cfg,hubnums,alldata_reduc );
plotusage( SATCOOR, CONNout, cfg, USAGE, DEPOTS, chhubs );


fprintf(' \n');




