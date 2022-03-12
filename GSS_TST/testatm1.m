function [ rho ] = testatm1( jd0 )
%ESTATM1 Summary of this function goes here
%   Detailed explanation goes here

secs = 86400*4;

jd0 = floor(jd0);
rho = zeros(secs,1);
jd = zeros(secs,1);

for t = 1:secs
    %fprintf('|');
    jd(t) = jd0+t/86400;
    [~,rtascsun,declsun] = astro.sun ( jd(t) ); %Vallado
    [month, day, ~] = astro.gdate(jd(t));
    [ F81,F107, Kp ] = astro.getsolar( jd(t), 3 );
    F107=F81; Kp = 3;
    thetag = astro.gstime(jd(t));
    [rho(t), rhoN, KK] = astro.atmgost2([6378.14+650 0 0], 650, (month-1)*30+day,thetag,rtascsun,declsun, F107, F81, Kp);
    
    if any(t==[157129 157130]) %&& 0
        fprintf('\n');
        fprintf('jd = %3.2f ; mm--dd = %3.0f--%3.2f ; rhoN = %g ; rho = %g \n',jd(t), month, day, rhoN, rho(t));
        fprintf('F107 = %3.2f; F81 = %3.2f, Kp = %3.2f \n',F107, F81, Kp);
        fprintf('rtascsun = %3.2f; declsun = %3.2f, theta_g = %3.2f \n',rtascsun,declsun, thetag);
        fprintf('KKs = %3.2f %3.2f %3.2f %3.2f %3.2f \n',KK(1),KK(2),KK(3),KK(4),KK(5));
    end
    
    
end

fprintf('\n');
x = 1:secs; x = x';%/3600;
figure; plot(x,rho); hold on; grid;


end

