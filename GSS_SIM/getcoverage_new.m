function [ gaguliki, fin, nsat, elevs_last, dists_last ] = getcoverage_new( alt, latgnd, longnd, elevcrit, plotmode, animate )
%GO_COVER Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    plotmode = 0;
    animate = false;
end

cdata = imread('EM_low.jpg');

gd0 = [2017 03 20 10 28 0]; % SET UP STARTTIME (UTC)

R0 = 6378.137;

jd0 = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6)); % find starttime Julian Date

mindthegap = 1; % if true, DO consider the gap
STEP = 90;

deg = pi/180;

a = R0+alt; e = 0; om = 0; in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
Om = 0*deg; u = 0*deg;

ini0 = [a e om Om in u];

t = 0:STEP:86400;
eph0 = astro.J2pert(ini0,t);
lt = length(t);

jt = jd0 + t/86400;

% gnd:

z0ecf = R0*sind(latgnd);
x0ecf = R0*cosd(latgnd)*cosd(longnd);
y0ecf = R0*cosd(latgnd)*sind(longnd);


plrange = 100; % INPUTS
sprange = 160;
plsize = plrange(end)-plrange(1)+1;
spsize = sprange(end)-sprange(1)+1;

gaguliki = zeros(plsize,spsize,lt);

finnan = ones(plsize,spsize);
nsat = ones(plsize,spsize);

alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));
beta = 90-elevcrit-alphacrit;

Rcov = R0 * sind(beta);
dcrit = R0*sind(beta)/sind(alphacrit); % critical swath distance

for i = plrange
    fprintf('Current planes = %g ->> ',i);
    im = int16(i-plrange(1) + 1);
    for j = sprange
        fprintf('%g > ',j);
        jm = int16(j-sprange(1) + 1);
        nsat(im,jm) = i*j;
        elevs_last = zeros(lt,i*j);
        dists_last = zeros(lt,i*j);
        [ ini0_cube, ~ , ~ ] = sat.seedconstel( ini0, jd0, i, j );
        [eph, ephcart] = astro.J2circ_cube(ini0_cube,t);
        for k = 1:lt
            %fprintf('|');
            thetag = astro.gstime(jt(k));
            
            %[ STATEN, CARTN, ~ ] = sat.seedconstel( eph0(k,:), jt(k), i, j );
            STATEN = permute(eph(:,k,:),[1,3,2]);
            CARTN = permute(ephcart(:,k,:),[1,3,2]);
            R0eci = math.R3(-thetag)*[x0ecf y0ecf z0ecf]';
            [MCOND, ~, elevs, dists] = sat.ifvisible(R0eci',CARTN(:,1:3),elevcrit);
            
            % bull^
            %bullshit = (1:i*j)';
            %bullshit = mod(bullshit(:),j) <= j/2; % drop half sats in plane
            %bullshit = mod((idivide(int32(bullshit(:)-1),j)+1),2) ~= 0; % drop planes
            bullshit = ones(i*j,1); % no bullshit
            MCOND = MCOND & bullshit;
            
            elevs_last(k,:) = elevs(:); elevs_last(k,~bullshit) = -90;
            dists_last(k,:) = dists(:); dists_last(k,~bullshit) = inf;
            
            
            if mindthegap  % we might be 'inthegap', it is then ignored
                LAN1 = STATEN(1,4) - thetag;
                LAN1 = mod(LAN1,2*pi);
                LDN1 = mod((LAN1+pi),2*pi);
                deltaAN = abs(LAN1/deg - longnd);
                deltaDN = abs(LDN1/deg - longnd);
                inthegapLAN = abs(deltaAN) < 20 || abs(deltaAN - 360) < 20;
                inthegapLDN = abs(deltaDN) < 20 || abs(deltaDN - 360) < 20;
                inthegap = inthegapLAN || inthegapLDN;
            else
                inthegap = false; % we are never 'inthegap', consider everything equally
            end
            
            
            
            gaguliki(im,jm, k) = sum(MCOND);
            if (j == 14 && i == 7) && k>=1 && animate
                %cla;
                vis_drawstate( CARTN, jt(k), i, j, cdata );
                add3d.all_coverage( CARTN(MCOND & (elevs< 15),1:3), alphacrit, ':w', 'c' );
                add3d.all_coverage( CARTN(MCOND & (elevs>=15),1:3), alphacrit, ':g', 'c' );
                add3d.majorcities(jt(k),1);
                add3d.sunlight(jt(k));
                add3d.moonpoint(jt(k));
                view(thetag/deg+longnd+90-15,latgnd+15);
                
                plot3(R0eci(1)*1e3,R0eci(2)*1e3,R0eci(3)*1e3,'dr', 'markerfacecolor', 'r'); % ground station
                plot3(CARTN(bullshit,1)*1e3,CARTN(bullshit,2)*1e3,CARTN(bullshit,3)*1e3,'or', 'markerfacecolor', 'c');
                plot3(CARTN(MCOND,1)*1e3,CARTN(MCOND,2)*1e3,CARTN(MCOND,3)*1e3,'or', 'markerfacecolor', 'g');
                %plot3(CARTN(1,1)*1e3,CARTN(1,2)*1e3,CARTN(1,3)*1e3,'*r', 'markerfacecolor', 'r');
                
                for l = 1:numel(MCOND)
                    if MCOND(l)
                        d = sqrt(sum( (R0eci'-CARTN(l,1:3)).^2,2 ));
                        hue = 0.5 - 0.5*(d-alt) / (dcrit-alt);
                        plot3([R0eci(1)*1e3 CARTN(l,1)*1e3],[R0eci(2)*1e3 CARTN(l,2)*1e3],[R0eci(3)*1e3 CARTN(l,3)*1e3],...
                            '-y', 'linewidth',2, 'color',hsv2rgb(hue,0.9,0.9));
                    end
                end
                drawnow;
                %pause(0.01);
            end
            
            
            
            %cond = k>2 && gaguliki(im,jm, k) < 1 && gaguliki(im,jm, k-1) < 1 && ~inthegap; % && gaguliki(im,jm, k-2) <= 1;
            cond = gaguliki(im,jm, k) < 1 && ~inthegap;
            %cond = 0; % ingore shit
            if cond
                finnan(im,jm) = nan;
                nsat(im,jm) = nan;
                break
            end
        end
        %CARTN(1,:)
    end
    fprintf('\n');
end

fin = mean(gaguliki,3) .* finnan;

%plot

switch plotmode
    case 1 % plot full visibility parameters
        x = (0:(lt-1))'*STEP / 60;
        
        elev_full = elevs_last;
        elev_full(elev_full < elevcrit) = nan;
        
        dist_full = dists_last;
        dist_full(dist_full > dcrit) = nan;
        
        figure;plot(x,elev_full);
        hold on;grid;xlabel('time, min');ylabel('elevation, deg');
        figure;plot(x,dist_full);
        hold on;grid;xlabel('time, min');ylabel('distance, km');
        plot(x,dcrit*ones(size(x)),'--r','linewidth',2);
        
    case 2 % plot optimal visibility parameters
        x = (0:(lt-1))'*STEP / 60;
        
        elev_desc = sort(elevs_last,2,'descend');
        elev_desc(elev_desc < elevcrit) = nan;
        
        dist_asc  = sort(dists_last,2,'ascend');
        dist_asc(dist_asc > dcrit) = nan;
        
        figure;plot(x,elev_desc(:,1));
        hold on;grid;xlabel('time, min');ylabel('elevation, deg');
        figure;plot(x,dist_asc(:,1), '-r');
        hold on;grid;xlabel('time, min');ylabel('distance, km');
        %plot(x,dcrit*ones(size(x)),'--r','linewidth',2);
        
        gagaga = elev_desc(:,1);
        figure; hold on;
        for q = 1:lt
            if gagaga(q) > 25
                hue = 0.5;
            elseif gagaga(q) > 10
                hue = 0.5*(gagaga(q)-10) / 15;
            else
                hue = 0;
            end
            plot([x(q) x(q)], [0 1], 'color',hsv2rgb([hue 0.9 0.9]), 'linewidth', 2);
            xlim ([135 355]);
        end

end


end

