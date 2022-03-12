function [ gaguliki, fin, nsat ] = GO_cover( alt, latgnd, longnd )
%GO_COVER Summary of this function goes here
%   Detailed explanation goes here

gd0 = [2017 03 20 10 28 0]; % SET UP STARTTIME (UTC)

%
% vtau = 8.192; % specify typeV session duration
% dtau = 0.192;  % specify typeN session duration
% stau = 2.048; % specify typeS session duration

R0 = 6378.137;

% byteperday1 = 76*512;
% byteperday2 = 450*12;
% byteperday3 = 226*128;

% br_m = 500;
% br_f = 500;

jd0 = jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6)); % find starttime Julian Date

mindthegap = 1;
STEP = 90;

deg = pi/180;

a = R0+alt; e = 0; om = 0; in = sunsyn(alt)*deg; % specify initial orbit parameters
Om = 0*deg; u = 0*deg;

ini0 = [a e om Om in u];
eph0 = J2pert(ini0,86400,STEP);
t = eph0(:,1); lt = length(t);
eph0 = eph0(:,2:7);

jt = jd0 + t/86400;

% gnd:

z0ecf = R0*sind(latgnd);
x0ecf = R0*cosd(latgnd)*cosd(longnd);
y0ecf = R0*cosd(latgnd)*sind(longnd);


plrange = 96:97; % INPUTS
sprange = 160:165;
plsize = plrange(end)-plrange(1)+1;
spsize = sprange(end)-sprange(1)+1;

gaguliki = zeros(plsize,spsize,lt);

finnan = ones(plsize,spsize);
nsat = ones(plsize,spsize);

elevcrit = 60; % deg
%alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));

for i = plrange
    fprintf('Current planes = %g ->> ',i);
    im = int16(i-plrange(1) + 1);
    for j = sprange
        fprintf('%g > ',j);
        jm = int16(j-sprange(1) + 1);
        nsat(im,jm) = i*j;
        for k = 1:lt
            %fprintf('|');
            [ STATEN, CARTN, ~ ] = seedconstel( eph0(k,:), jt(k), i, j );
            R0eci = R3(-gstime(jt(k)))*[x0ecf y0ecf z0ecf]';
            MCOND = ifvisible(R0eci',CARTN(:,1:3),elevcrit);
            %             R0ECIN = repmat(R0eci',i*j,1);
            %             R01N = CARTN(:,1:3) - R0ECIN;
            %             dists = sqrt(sum(R01N.*R01N,2));
            %             rads = sqrt(sum(CARTN(:,1:3).*CARTN(:,1:3),2));
            %             DOTN = dot(R01N,CARTN(:,1:3),2);
            %             cosalphas = DOTN ./ dists ./rads;
            %             alphas = acosd(cosalphas);
            %             MCOND = ( alphas < alphacrit & dot(R01N,R0ECIN,2) > 0 );
            
            
            
            if mindthegap
                LAN1 = STATEN(1,4) - gstime(jt(k));
                LAN1 = mod(LAN1,2*pi);
                LDN1 = mod((LAN1+pi),2*pi);
                deltaAN = abs(LAN1/deg - longnd);
                deltaDN = abs(LDN1/deg - longnd);
                inthegapLAN = abs(deltaAN) < 20 || abs(deltaAN - 360) < 20;
                inthegapLDN = abs(deltaDN) < 20 || abs(deltaDN - 360) < 20;
                inthegap = inthegapLAN || inthegapLDN;
            else
                inthegap = false;
            end
            
            
            
            gaguliki(im,jm, k) = sum(MCOND);
            if (j == 35 && i == 15) && k==1095
                vis_drawstate( CARTN, jt(k), i, j );
                plot3(R0eci(1)*1e3,R0eci(2)*1e3,R0eci(3)*1e3,'dr', 'markerfacecolor', 'r');
                plot3(CARTN(MCOND,1)*1e3,CARTN(MCOND,2)*1e3,CARTN(MCOND,3)*1e3,'or', 'markerfacecolor', 'g');
                plot3(CARTN(1,1)*1e3,CARTN(1,2)*1e3,CARTN(1,3)*1e3,'*r', 'markerfacecolor', 'r');
            end
            
            
            
            %cond = k>2 && gaguliki(im,jm, k) <= 1 && gaguliki(im,jm, k-1) <= 1 ; % && gaguliki(im,jm, k-2) <= 1;
            cond = gaguliki(im,jm, k) < 1 && ~inthegap;
            %cond = 0; % ingore shit
            if cond
                finnan(im,jm) = nan;
                nsat(im,jm) = nan;
                break
            end
        end
    end
    fprintf('\n');
end

fin = mean(gaguliki,3) .* finnan;

end

