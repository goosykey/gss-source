function [ Nreal, Ntot, Ndes, Nass, POWER ] = CYCLO_SIM( EPHEMERIS, netdata, jd, map, mapdata, BIGASSMASK  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fprintf(char(datetime('now')));
fprintf(' > Calculation begins. \n');

t = EPHEMERIS(:,1);
NN = numel(t);

Nreal = zeros(NN,3);
Ndes = zeros(NN,3);
Nass = zeros(NN,3);
Ntot = zeros(NN,3);
POWER = zeros(NN,1);

vtau = netdata{4};
dtau = netdata{5};
stau = netdata{6};

jt = jd + t/86400;

period = t(end);
tick = t(2)-t(1);
hours = ceil(period/3600);

[ SESDATA, SESCOORD, SESBEG ] = gettrafficearth( jd, 1800, map, mapdata, BIGASSMASK, netdata );

global temp4 temp5 temp6 temp7;

for i = 1:NN
    ti = jt(i);
    if max(SESBEG) < ti
        fprintf(char(datetime('now')));
        fprintf(' > Calculation suspended. Progress: %g%% \n', (ti-jd)/(jt(end)-jd)*100);
        fprintf('Session data not enough. Generating new session data...\n ');
        [ sesdata, sescoord, sesbeg ] = gettrafficearth( ti, 1800, map, mapdata, BIGASSMASK, netdata );
        BIGSHIT = [double([SESDATA;sesdata]),double([SESCOORD;sescoord]),[SESBEG;sesbeg]]; % in order to sort by session type !IMPORTANT
        BIGSHIT = sortrows(BIGSHIT,1);
        
        SESDATA = int8(BIGSHIT(:,1:2));
        SESCOORD = BIGSHIT(:,3:4);
        SESBEG = BIGSHIT(:,5);
        fprintf('Done. Calculation resumed.\n');
    end
    
    fprintf(' -> Progress: %3.2f%%', (ti-jd)/(jt(end)-jd)*100);
    
    INTOSESSION = round((ti-SESBEG)*86400);
    SESBEG((INTOSESSION(:) > vtau) & (SESDATA(:,1)==1),:) = nan; % v.sessions already over
    SESBEG((INTOSESSION(:) > dtau) & (SESDATA(:,1)==2),:) = nan; % d.sessions already over
    SESBEG((INTOSESSION(:) > stau) & (SESDATA(:,1)==3)) = nan; % s.sessions already over
    
    SESDATA(isnan(SESBEG),:) = [];
    SESCOORD(isnan(SESBEG),:) = [];
    SESBEG(isnan(SESBEG)) = [];
    
    IS = round((ti-SESBEG)*86400);
    SD = SESDATA(IS>=0,:);
    SC = SESCOORD(IS>=0,:);
    SB = SESBEG(IS>=0);
    
    if i == 106
        temp4 = IS;
        temp5 = SD;
        temp6 = SC;
        temp7 = SB;
    end
    
    
    a = EPHEMERIS(i,2);
    e = EPHEMERIS(i,3);
    om = EPHEMERIS(i,4);
    Om = EPHEMERIS(i,5);
    in = EPHEMERIS(i,6);
    u = EPHEMERIS(i,7);
    
    %numel(SESBEG)
       
    [ Nreali, Pout, Ntoti, Ndesi, Nassi] = getener( a, e, om, Om, in, u, netdata, SD, SC, ti );
    fprintf(' -> Pout = %6.2f W -> Zone : %04.0f -> Distr : %04.0f -> Add.assist : %04.0f -> Denied : %04.0f \n',...
        Pout, sum(Ntoti), sum(Ndesi), sum(Nassi)-sum(Nreali), sum(Ndesi)-sum(Nassi) );
    
    
    Nreal(i,:) = Nreali;
    Ndes(i,:) = Ndesi;
    Nass(i,:) = Nassi;
    Ntot(i,:) = Ntoti;
    POWER(i) = Pout;
end

end

