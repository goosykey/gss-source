function [ output_args ] = supermegaplot(t, y, yfil, opts)
%SUPERMEGAPLOT
% use for visualization
%
%   by A.Kharlan, © Yaliny, 2015
%
%   Options:
%       'a' - plot semi-major axis
%       'e' - plot eccentricity
%       'o' - plot arg. of perigee
%       'O' - plot RAAN
%       'i' - plot inclonation
%       'u' - plot arg.of latitude
%       'v' - plot arg. of latitude -> [0..2pi)
%       'm' - plot mass
%       'r' - plot Earth magnitude
%       'b' - plot sun Beta angle
%       'n' - plot true anomaly
%       'p' - plot focal parameter
%       'h' - plot WGS84 altitude
%       'V' - show mean value on the whole interval
%       'F' - fullscreen charts
% 
% 
%     INPUTS:
%     t, y      : "true" ephemeris
%     t, yfil   : "mean" ephemeris
%     opts      : options (see above)

%Arg management
if nargin < 4
    opts = 'aeimu';
end

if nargin < 3
    plotyavg = false;
else
    plotyavg = true;
end

global STARTTIME;

optV = ~isempty(strfind(opts,'V'));
optF = ~isempty(strfind(opts,'F'));

if yfil == 0
    yfil = y;
end

p = y(:,1);
pfil = yfil(:,1);
l1 = y(:,2);
l1fil = yfil(:,2);
l2 = y(:,3);
l2fil = yfil(:,3);
e = sqrt(l1.^2 + l2.^2);
efil = sqrt(l1fil.^2 + l2fil.^2);
a = p./(1-e.^2);
afil = pfil./(1-efil.^2);
om = perigee(l1,l2);
omfil = perigee(l1fil,l2fil);
Om = y(:,4);
Omfil = yfil(:,4);
in = y(:,5);
infil = yfil(:,5);
u = y(:,6);
v = u;

if ~isempty(strfind(opts,'v')) || ~isempty(strfind(opts,'r')) || ~isempty(strfind(opts,'h'))
    v = mod(v,2*pi);
%     while sum(v(:)>2*pi) ~= 0
%         v(v(:)>=2*pi) = v(v(:)>=2*pi) - 2*pi;
%     end
end

%om = perigee(l1,l2);

au1 = 149597871; % 1 Astronomic Unit in kms
mu = 398600.44;

t = t ./ 86400;

bet = zeros(length(t),1);

scrsz = get(groot,'ScreenSize');
if optF
    scrH = scrsz(4)-30;
    scrW = scrsz(3);
    scrMx = [1,30; 1,30; 1,30; 1,30];
else
    scrH = scrsz(4)/2;
    scrW = scrsz(3)/2;
    scrMx = [1,scrH; scrW,scrH; 1,30; scrW,30];
end

% SEMI-MAJOR
if ~isempty(strfind(opts,'a'))
    figure('Name','Report: Semi-major axis','NumberTitle','off', 'OuterPosition',[scrMx(1,:), scrW, scrH]);
    plot(t, a, 'r', 'LineWidth', 2); grid; hold; % Построение графика и сетки
    plot(t, afil, 'y', 'LineWidth', 2);
    if optV
        pa = polyfit(t,a,1);
        aavg = polyval(pa,t);
        plot(t, aavg, 'c', 'LineWidth', 2);
        legend('a(t)','a_f_i_l(t)','a_a_v_g(t)');
    else 
        legend('a(t)','a_f_i_l(t)');
    end
    xlabel('time, days');
    ylabel('semi-major axis, km');
end


% INCLINATION
if ~isempty(strfind(opts,'i'))
    figure('Name','Report: inclination','NumberTitle','off', 'OuterPosition',[scrMx(2,:), scrW, scrH]);
    plot(t, in .*180 ./pi, 'LineWidth', 2); grid; hold;
    plot(t, infil.*180 ./pi, 'y', 'LineWidth', 2);
    if optV
        pin = polyfit(t,in,1);
        inavg = polyval(pin,t);
        plot(t, inavg .*180 ./pi, 'c', 'LineWidth', 2);
        legend('i(t)','i_f_i_l(t)','i_a_v_g(t)');
    else 
        legend('i(t)','i_f_i_l(t)');
    end
    xlabel('time, days');
    ylabel('inclination, deg');
end

% RAAN
if ~isempty(strfind(opts,'O'))
    figure('Name','Report: RAAN','NumberTitle','off');
    plot(t, Om*180/pi, 'LineWidth', 2); grid; 
    legend('O(t)');
    xlabel('time, days');
    ylabel('RAAN, deg');
end

% ECCENTRICITY
if ~isempty(strfind(opts,'e'))
    figure('Name','Report: Eccentricity','NumberTitle','off', 'OuterPosition',[scrMx(3,:), scrW, scrH]);
    plot(t, e, 'm', 'LineWidth', 2); grid; hold;
    plot(t, efil, 'y', 'LineWidth', 2);
    legend('e(t)');
    if optV
        pe = polyfit(t,e,1);
        eavg = polyval(pe,t);
        plot(t, eavg, 'g', 'LineWidth', 2);
        legend('e(t)','e_f_i_l(t)','e_a_v_g(t)');
    else
        legend('e(t)','e_f_i_l(t)');
    end
    xlabel('time, days');
    ylabel('eccentricity');
end

% ARG OF PERIGEE
if ~isempty(strfind(opts,'o'))
    figure('Name','Report: Arg. of perigee','NumberTitle','off', 'OuterPosition',[scrMx(4,:), scrW, scrH]);
    plot(t, om, 'c', 'LineWidth', 2); grid; hold;
    plot(t, omfil, 'r', 'LineWidth', 2);
    legend('o(t)');
    if optV
        pom = polyfit(t,om,1);
        omavg = pom(1)*t + pom(2);
        plot(t, omavg, 'g', 'LineWidth', 2);
        legend('o(t)','o_f_i_l(t)','o_a_v_g(t)');
    else
        legend('o(t)','o_f_i_l(t)');
    end
    xlabel('time, days');
    ylabel('arg. of perigee, rad');
end

% ARG OF LATITUDE
if ~isempty(strfind(opts,'u'))
    figure('Name','Report: Arg. of latitude','NumberTitle','off');
    plot(t, u, 'r', 'LineWidth', 2); grid; hold;
    T0 = 2*pi*sqrt(a(1)^3/mu);
    unom = u(1) + 2*pi*(t*86400 - t(1))/T0;
    plot(t, unom, '--c', 'LineWidth', 2);
    legend('u(t)', 'u_n_o_m(t)');
    xlabel('time, days');
    ylabel('arg. of latitude, rad');
end

% ARG OF LATITUDE
if ~isempty(strfind(opts,'l'))
    figure('Name','Report: Lambdas','NumberTitle','off');
    plot(t, l1, 'b', 'LineWidth', 2); grid; hold;
    plot(t, l1fil, 'y', 'LineWidth', 2);
    plot(t, l2, 'g', 'LineWidth', 2);
    plot(t, l2fil, 'm', 'LineWidth', 2);
    legend('\lambda_1(t)', '\lambda_1_f_i_l(t)', '\lambda_2(t)', '\lambda_2_f_i_l(t)');
    xlabel('time, days');
    ylabel('\lambda_1, \lambda_2');
end

% ARG OF LAT DIMINISHED
if ~isempty(strfind(opts,'v'))
    figure('Name','Report: Arg. of latitude','NumberTitle','off');
    plot(t, v, 'LineWidth', 2); grid; 
    legend('u(t)');
    xlabel('time, days');
    ylabel('arg. of latitude, rad');
end

% MASS
if ~isempty(strfind(opts,'m')) 
    figure('Name','Report: MASS','NumberTitle','off');
    plot(t, y(:,7), 'LineWidth', 2); grid;
    legend('M(t)');
    xlabel('time, days');
    ylabel('Mass, kg');
end

% SUN BETA ANGLE
if ~isempty(strfind(opts,'b'))
    figure('Name','Report: Beta-angle','NumberTitle','off');
    ST0 = STARTTIME;
    jd0 = jday(ST0(1), ST0(2), ST0(3), ST0(4), ST0(5), ST0(6));
    jd = t + jd0;
    for q = 1:length(t)
        [rsun,rtascsun,declsun] = sun ( jd(q) );
        %bet(q) = asin(cos(declsun)*sin(y(q,5))*sin(y(q,4)-rtascsun)+sin(declsun)*cos(y(q,5)));
        normxyz = tsw2xyz([0;0;-1],Om(q),in(q),u(q));
        sunxyz = rsun.*au1;
        fif = acos(sunxyz*normxyz/(norm(sunxyz)*norm(normxyz)));
        bet(q) = pi/2 - fif;
        
    end
    plot(t, bet.*180./pi, 'LineWidth', 2); grid;
    legend('BETA(t)');
    xlabel('time, days');
    ylabel('Beta-angle, deg');
end

if ~isempty(strfind(opts,'r')) || ~isempty(strfind(opts,'h'))
    figure('Name','Report: Radius','NumberTitle','off');
    %nu = v-om;
    %nu(nu(:)<0) = nu(nu(:)<0) + 2*pi;
    %r = p .* (1 + e.*cos(nu));
    [R1,V1] = randv(a,e,in,Om,zeros(length(v),1),v);
    R=R1'; V=V1';
    r = sqrt(R(:,1).^2+R(:,2).^2+R(:,3).^2);
    if ~isempty(strfind(opts,'r'))
        plot(t, r, 'LineWidth', 2); grid;
        xlabel('time, days');
        ylabel('r, km');
    end
end

if ~isempty(strfind(opts,'h'))
    figure('Name','Report: altitude (WGS84)','NumberTitle','off');
    nu = v-om;
    nu(nu(:)<0) = nu(nu(:)<0) + 2*pi;
    r = p .* (1 + e.*cos(nu));
    L = length(r);
    rxyz = R;
    hvec = zeros (L,3);
    h = zeros (L,1);
    rdot =  zeros(L,3);
    awgs84 = 6378.137;
    bwgs84 = 6356.752;
    for jj = 1:L
        rdot(jj,2) = sqrt(1/(rxyz(jj,1)^2/(awgs84^2*rxyz(jj,2)^2)+1/awgs84^2+rxyz(jj,3)^2/(bwgs84^2*rxyz(jj,2)^2)));
        rdot(jj,1) = rdot(jj,2)*rxyz(jj,1)/rxyz(jj,2);
        rdot(jj,3) = rdot(jj,2)*rxyz(jj,3)/rxyz(jj,2);
        %norm(rdot(jj,:))
        %hvec(jj,:) = rxyz(jj,:)-rdot(jj,:);
        %h(jj) = norm(hvec(jj,:));
        h(jj) = norm(rxyz(jj,:)) - norm(rdot(jj,:));
    end
    plot(t, h, 'LineWidth', 2); grid;
    xlabel('time, days');
    ylabel('h, km');
end


%clc;

end

