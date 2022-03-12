function [ t,y,yavg ] = GOFULL( t_fin, dayperiod, h, Y, opts)
%PRIMARY EPHEMERIS GENERATOR
%
%   by A.Kharlan, © Yaliny, 2015
%
%   Options:
%       'P' - use 'Precise' motion model:
%           J20...J44+prec&nut, Sun, Moon, Atmosphere;
%       'A' - use 'Advanced' motion model:
%           J2, J4, Sun, Moon, Atmosphere;
%       'S' - use 'Simplified' motion model (DEFAULT):
%           J2, Atmosphere;
%

tic;

if nargin <5 
    opts = '';
end

optP = ~isempty(strfind(opts,'P'));
optA = ~isempty(strfind(opts,'A'));
if optP
    prec = '_precise';
elseif optA
    prec = '_advanced';
else
    prec = '';
end

t = zeros(10000000,1);
y = zeros(10000000,7);
yavg = zeros(10000000,7);

per = dayperiod * 86400;
t0 = 0;
icur = 1;
skmade = 0;
y0 = Y;
options = odeset('MaxStep', 240);
anom = y0(1);
%anom = 7178.142;
ak = y0(1);
da = 0.5;
%da = 500;

mu = 398600.44;

scrsz = get(groot,'ScreenSize');
figure('Name','Report: FOCAL','NumberTitle','off', 'OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
grid; hold;
xlabel('time, days');
ylabel('Focal parameter, km');

subint = ['earth353',prec];
TOCtotal = 0;

while t0 < t_fin
    daypassed = t0/86400;
    clc;
    fprintf('day = %g ; SK sessions = %g \n', daypassed, skmade);
    TOC=toc;
    TOCtotal = TOCtotal+TOC;
    fprintf('last step = %4.2f sec; time elapsed = %4.2f min\n', TOC, TOCtotal/60);
    fprintf('%3.2f%% COMPLETE;', t0/t_fin*100);
    tic;
    if anom - ak < da
        clear ti yi;
        [ti, yi] = ode45(subint, [t0:h:t0+per], y0, options);
        Lcur = length (ti);
        plot(ti/86400, yi(:,1), 'r', 'LineWidth', 2);
        t(icur:icur+Lcur-1) = ti;
        y(icur:icur+Lcur-1,:) = yi;
        pavg = zeros(1,7);
        for ii = 1:7
            pavg(ii) = (max(yi(:,ii))+min(yi(:,ii)))/2;
        end
        yia = yi;
        for jj = 1:length(ti)
            yia(jj,:) = pavg;
        end
        plot(ti/86400, yia(:,1), 'c', 'LineWidth', 2);
        yavg(icur:icur+Lcur-1,:) = yia;
        pk = pavg(1);
        l1k = pavg(2);
        l2k = pavg(3);
        ek = sqrt(l1k^2 + l2k^2);
        ak = pk / (1 - ek^2);
               
        t0 = t0 + per;
        y0 = yi(length(yi(:,1)),:);
        icur = icur+Lcur;
    else
        skmade = skmade+1;
        while ak-anom<da
            clear ti yi;
            T = 2*pi * sqrt(y0(1)^3/mu);
            [ti, yi] = ode45([subint,'_Tprop'], [t0:h:t0+T], y0, options);
            Lcur = length (ti);
            plot(ti/86400, yi(:,1), 'b', 'LineWidth', 2);
            t(icur:icur+Lcur-1) = ti;
            y(icur:icur+Lcur-1,:) = yi;
            pavg = zeros(1,7);
            for ii = 1:7
                pavg(ii) = (max(yi(:,ii))+min(yi(:,ii)))/2;
            end
            yia = yi;
            for jj = 1:length(ti)
                yia(jj,:) = pavg;
            end
            plot(ti/86400, yia(:,1), 'c', 'LineWidth', 2);
            yavg(icur:icur+Lcur-1,:) = yia;
            pk = pavg(1);
            l1k = pavg(2);
            l2k = pavg(3);
            ek = sqrt(l1k^2 + l2k^2);
            ak = pk / (1 - ek^2);
        
            t0 = t0 + T;
            y0 = yi(length(yi(:,1)),:);
            icur = icur+Lcur;
        end %endwhile
        
    end %endif anom - ak < da
    
end %endwhile t0 < t_fin

t = t(y(:,1)~=0);
y = y(y(:,1)~=0,:);
yavg = yavg(y(:,1)~=0,:);

end

