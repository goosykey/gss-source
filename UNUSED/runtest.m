%=== ISL Traffic load testing

%   Specify [ <Number of orbits> , <Number of satellites per orbit> ]
cfg = [8,21];

%   Failure probability
prob = 0.1;

%   Generate the ISL net with flawed connections
[ SATCOOR, CONNout ] = generatenet( cfg, prob);

%   ANALYZE LOAD. Arg. #3 = number of gnd stations. Also, optionally
%   include depots manual input:
%   [ DEPOTS, USAGE, isol ] = ISCLoad( CONNout,cfg,GNDSTAT, [depot1, depot2, depot3...] );
[ DEPOTS, USAGE, isol, chhubs ] = ISCLoad( CONNout,cfg,40,[],9308,0.3 );

%   plot results
plotusage( SATCOOR, CONNout, cfg, USAGE, DEPOTS, chhubs );