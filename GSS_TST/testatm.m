function [ rho, F81 ] = testatm( jd0 )
%TESTATM Summary of this function goes here
%   Detailed explanation goes here

jd0 = round(jd0)+0;

td = 365*15;
rho = zeros(td,1);
rhoN = zeros(td,1);
KK = zeros(td,5);
jd = zeros(td,1);
F81 = zeros(td,1);
[month0, day0, year0] = astro.gdate(jd0);

for q = 1:td
    %fprintf('|');
    jd(q) = jd0+q;
    [~,rtascsun,declsun] = astro.sun ( jd0+q ); %Vallado
    [month, day, ~] = astro.gdate(jd0+q);
    [ F81(q),F107, Kp ] = astro.getsolar( jd0+q, 3 );
    thetag = astro.gstime(jd0+q);
    [rho(q), rhoN(q), KK(q,:)] = astro.atmgost2([6378.14+650 0 0], 650, (month-1)*30+day,thetag,rtascsun,declsun, F107, F81(q), Kp);
    
    if any(round(jd(q))==[2460294 2460295 2464346 2464347]) %&& 0
        fprintf('\n');
        fprintf('jd = %3.2f ; mm--dd = %3.0f--%3.2f ; rhoN = %g ; rho = %g \n',jd(q), month, day, rhoN(q), rho(q));
        fprintf('F107 = %3.2f; F81 = %3.2f, Kp = %3.2f \n',F107, F81(q), Kp);
        fprintf('rtascsun = %3.2f; declsun = %3.2f, theta_g = %3.2f \n',rtascsun,declsun, thetag);
        fprintf('KKs = %3.2f %3.2f %3.2f %3.2f %3.2f \n',KK(q,1),KK(q,2),KK(q,3),KK(q,4),KK(q,5));
    end

end

fprintf('\n');
x = 1:td; x = x'/365 + year0 + (month0-1)/12 + day0/365;
%x = (1:td)';
figure; plot(x,rho); hold on; grid;
%figure; plot(x,F81); hold on; grid;

end

