function [ t,y,yfil ] = PEG( t_0, h, t_fin, dayperiod, Y0, Yn, opts)
%PRIMARY EPHEMERIS GENERATOR
%
%   by A.Kharlan, © Yaliny, 2015
%
%   Options:
%       'P' - use 'Precise' motion model:
%           J20...J88+prec&nut, Sun, Moon, Atmosphere;
%       'A' - use 'Advanced' motion model:
%           J2, J4, Sun, Moon, Atmosphere;
%       'S' - use 'Simplified' motion model (DEFAULT):
%           J2, Atmosphere;
%
%     INPUTS:
%     t_0       : start time, sec
%     h         : ephemeris step, sec
%     t_fin     : stop time, sec
%     dayperiod : check period, days
%     Y0        : starting point (parameters at t_0), [km,1,1,rad,rad,rad,kg]
%     Yn        : nominal parameters, [km,1,1,rad,rad,rad,kg]
%     opts      : options (see above)
%
%     OUTPUTS:
%     t, y      : "true" ephemeris
%     t, yfil   : "mean" ephemeris
%
%     GLOBAL VARIABLES USED:
%     PROPXSI  : yaw XSI-angle ports to the motion model
%     PREC     : defined by 'OPTIONS', this determines the modelling precision
%     STAGESTR : current calculation stage description
%     BOUND    : defined here, this is sin(ulim) to limit active
%               stage during each orbit. Negative value means phasing,
%               positive value means stationkeeping;

tic;

global PREC PROPXSI STAGESTR BOUND SATNAME;

if nargin <5
    opts = '';
end

optP = ~isempty(strfind(opts,'P'));
optA = ~isempty(strfind(opts,'A'));
if optP
    PREC = 2;
elseif optA
    PREC = 1;
else
    PREC = 0;
end

t = zeros(10000000,1);
y = zeros(10000000,7);
yfil = zeros(10000000,7);

per = dayperiod * 86400;
t0 = t_0;
icur = 1;
skmade = 0;
y0 = Y0;
yn = Yn;
options = odeset('MaxStep', 240,'RelTol', 1e-05);
anom = yn(1);
%anom = 7178.142;
ak = y0(1);

if strcmp (SATNAME, 'beacon')
    da = 10000;
else
    da = 0.5;
end
%da = 500;

mu = 398600.44;
col = 0;

TOCtotal = 0;

while t0 < t_fin
    daypassed = t0/86400;
    clc;
    
    if ak < 6500
        col = 1;
        break
    end
    
    fprintf(STAGESTR);
    fprintf('day = %g ; SK sessions = %g \n', daypassed, skmade);
    TOC=toc;
    TOCtotal = TOCtotal+TOC;
    fprintf('last step = %4.2f sec; time elapsed = %4.2f min\n', TOC, TOCtotal/60);
    fprintf('%3.2f%% COMPLETE;\n', (t0-t_0)/(t_fin-t_0)*100);
    tic;
    if anom - ak < da
        clear ti yi;
        PROPXSI = -1;
        [ti, yi] = ode45('earth353_FULL', [t0:h:t0+per], y0, options);
        Lcur = length (ti);
        %plot(ti/86400, yi(:,1), 'r', 'LineWidth', 2);
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
        %plot(ti/86400, yia(:,1), 'c', 'LineWidth', 2);
        yfil(icur:icur+Lcur-1,:) = yia;
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
            PROPXSI = 0;
            BOUND = 0.984808;
            BOUND = 0;
            [ti, yi] = ode45('earth353_FULL', [t0:h:t0+T], y0, options);
            Lcur = length (ti);
            %plot(ti/86400, yi(:,1), 'b', 'LineWidth', 2);
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
            %plot(ti/86400, yia(:,1), 'c', 'LineWidth', 2);
            yfil(icur:icur+Lcur-1,:) = yia;
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
yfil = yfil(y(:,1)~=0,:);

t = t(t(:)<=t_fin);
y = y(t(:)<=t_fin,:);
yfil = yfil(t(:)<=t_fin,:);

daypassed = t0/86400;
clc;
fprintf(STAGESTR);
fprintf('day = %g ; SK sessions = %g \n', daypassed, skmade);
TOC=toc;
TOCtotal = TOCtotal+TOC;
fprintf('last step = %4.2f sec; time elapsed = %4.2f min\n', TOC, TOCtotal/60);
if col
    fprintf('WARNING: the satellite has collapsed.\n');
else
    fprintf('%3.2f%% COMPLETE;\n', (t0-t_0)/(t_fin-t_0)*100);
end

end

