function [ t,y,yfil ] = PEG( jd0, h, delta_t, dayperiod, Y0, Yn, opts)
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

global SATDATA PREC STAGESTR;

if nargin <5
    opts = '';
end

optP = contains(opts,'P');
optA = contains(opts,'A');
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
t0 = 0;
icur = 1;
skmade = 0;
y0 = Y0;
yn = Yn;
options = odeset('MaxStep', 240,'RelTol', 1e-05);
anom = yn(1);
ak = y0(1);


da = 1000; % CORREC
BOUNDDEG = 45;


mu = 398600.44;
col = 0;

TOCtotal = 0;

while t0 < delta_t
    daypassed = t0/86400;
    clc;
    
    if ak < 6500
        col = 1;
        break
    end
    
    fprintf(STAGESTR);
    if anom - ak < da
        fprintf('day = %g ; SK sessions done = %g \n', daypassed, skmade);
    else
        fprintf('day = %g ; STATION KEEPING #%g NOW! \n', daypassed, skmade+1);
    end
    TOC=toc;
    TOCtotal = TOCtotal+TOC;
    fprintf('last step = %4.2f sec; time elapsed = %4.2f min\n', TOC, TOCtotal/60);
    fprintf('%3.2f%% COMPLETE;\n', t0/delta_t*100);
    tic;
    if anom - ak < da
        %ak-6378.14
        clear ti yi;
        SATDATA{3} = zeros(360,1); % XSI PROP
        SATDATA{4} = zeros(360,1); % PROP FLAGS
        [ti, yi] = ode45('earth353_FULL', (t0+jd0*86400):h:(t0+per+jd0*86400), y0, options);
        ti = ti-jd0*86400; % jday to t
        Lcur = length (ti);
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
            %ak-6378.14
            clear ti yi;
            T = 2*pi * sqrt(y0(1)^3/mu);
            SATDATA{3} = zeros(360,1); % XSI PROP
            SATDATA{4} = phasing.getflgmass(BOUNDDEG,1);
            [ti, yi] = ode45('earth353_FULL', (t0+jd0*86400):h:(t0+T+jd0*86400), y0, options);
            ti = ti-jd0*86400; % jday to t
            Lcur = length (ti);
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

t = t(t(:)<=delta_t);
y = y(t(:)<=delta_t,:);
yfil = yfil(t(:)<=delta_t,:);

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
    fprintf('%3.2f%% COMPLETE;\n', t0/delta_t*100);
end

end

