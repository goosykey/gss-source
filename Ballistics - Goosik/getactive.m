function [ XSIact, deltaOmega, tact, yact, yfilact ] = getactive( Yn, Da, Di )
%ENGINE BURN STAGE corrector 
%   active phasing calculator
%   RAAN INDETERMINACY IS IGNORED IN THIS STAGE
%
%   by A.Kharlan, © Yaliny, 2015
%
%     INPUTS:
%     Yn : nominal orbit parameters  (BEGINS WITHOUT delta in RAAN)
%     Da : phasing semi-major axis difference
%     Di : phasing inclination difference
% 
%     OUTPUTS:
%     XSIact          : yaw XSI-angle for attitude thrust control
%     deltaOmega      : RAAN drift over the active stage
%     tact, yact      : "true" ephemeris for active stage
%     tact, yfilact   : "mean" ephemeris for active stage
% 
%     GLOBAL VARIABLES USED:
%     PROPXSI : yaw XSI-angle defined here ports to the motion model
%     PREC    : defined elsewhere, this determines the modelling precision
%     BOUND     : defined elsewhere, this is sin(ulim) to limit active
%               stage during each orbit. Negative value means phasing,
%               positive value means stationkeeping;

tic;
global PROPXSI BOUND STAGESTR;
STAGESTR = 'CALCULATING ACTIVE STAGE PARAMETERS...\n';
BOUND = -0.819152;

mu = 398600.44;
deg = pi/180;

p = Yn(1);
l1 = Yn(2);
l2 = Yn(3);
Om = Yn(4);
in = Yn(5);
u = Yn(6);
m = Yn(7);

en = sqrt(l1^2 + l2^2);

an = p * (1-en^2);
Dp = Da * (1-en^2);

orb = [an-Da en 0 10.9*pi/180 in-Di*deg pi/4 60];
y00 = orbitize(orb);


options = odeset('MaxStep', 240);
h = 10;

%ITERATION PARAMETERS (INIT)
itol = 0.005*deg;
stepno = 1;
maxstep = 100;
converged = false;
curstep = 1*deg;
goinup = true;
XSI = 80 * deg;
di = NaN;
t0 = NaN;

TOCtotal = 0;

while ~converged && stepno < maxstep;
    
    clc;
    fprintf('Da = %g ; Di = %g \n', Da, Di);
    TOC=toc;
    TOCtotal = TOCtotal+TOC;
    
    fprintf(STAGESTR);
    fprintf('last step = %4.2f sec; time elapsed = %4.2f min\n', TOC, TOCtotal/60);
    fprintf('Step # %g of %g;  = %4.2f deg\n', stepno, maxstep, curstep / deg);
    fprintf('current XSI = %4.2f deg \n', XSI / deg);
    fprintf('--------------------\n');
    fprintf('last counted duration %4.2f days; \n', t0/86400);
    fprintf('last counted residual "di" = %4.3f deg \n', di / deg);
    %fprintf('ink = %4.2f; in = %4.2f \n', ink/deg, in/deg);
    tic;
    
    PROPXSI = XSI;
    y0 = y00;
    t0 = 0;
    icur = 1;
    ak = y0(1);
    
    clear t y yavg;
    t = zeros(10000000,1);
    y = zeros(10000000,7);
    yavg = zeros(10000000,7);
   
    while ak < an
        clear ti yi;
        T = 2*pi * sqrt(y0(1)^3/mu);
        [ti, yi] = ode45('earth353_FULL', [t0:h:t0+T], y0, options);
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
    
    ink = pavg(5);
    di = ink - in;
    if abs(di) < itol
        converged = true;
    else
        if (di < 0 && ~goinup) || (di > 0 &&goinup)
            curstep = -curstep/2;
            goinup  = ~goinup;
        end               
    end
    
    stepno = stepno + 1;
    XSI = XSI + curstep;    
        
end

stepno = stepno - 1;
XSI = XSI - curstep;

fprintf('--------------------\n');
if converged
    fprintf('*CONVERGED in %g steps; XSIact = %4.2f deg\n', stepno, XSI / deg);
else
    fprintf('*NOT CONVERGED in %g steps; closest XSIact = %4.2f deg\n', stepno, XSI / deg);
end


fprintf('*last counted duration %4.2f days; \n', t0/86400);
fprintf('*last counted residual "di" = %4.3f deg \n', di / deg);

t = t(y(:,1)~=0);
y = y(y(:,1)~=0,:);
yavg = yavg(y(:,1)~=0,:);

% figure('Name','Report: FOCAL, LAST STEP','NumberTitle','off');
% grid; hold on;
% plot(t/86400, y(:,1), 'r', 'LineWidth', 2);
% plot(t/86400, yavg(:,1), 'c', 'LineWidth', 2);
% 
% figure('Name','Report: INCLIN, LAST STEP', 'NumberTitle','off');
% grid; hold on;
% plot(t/86400, y(:,5)/deg, 'b', 'LineWidth', 2);
% plot(t/86400, yavg(:,5)/deg, 'y', 'LineWidth', 2);c

tact = t;
yact = y;
yfilact = yavg;
XSIact = XSI;
deltaOmega = yavg (length(yavg),4)-yavg (1,4);

end

