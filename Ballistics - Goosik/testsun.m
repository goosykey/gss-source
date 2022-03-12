function [ MAM ] = testsun( h, Y, options )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global STARTTIME;
au1 = 149597871; % 1 Astronomic Unit in kms
deg = pi/180;

[ti, yi] = ode45('earth353_FULL', [0:h:86400], Y, options);

p = yi(:,1);
l1 = yi(:,2);
l2 = yi(:,3);
e = sqrt(l1.^2 + l2.^2);
a = p./(1-e.^2);
om = perigee(l1,l2); om(1)=0;
Om = yi(:,4);
in = yi(:,5);
u = yi(:,6);

[R,V] = randv(a,e,in,Om,om,u-om);

jd = jday(STARTTIME(1),STARTTIME(2),STARTTIME(3),STARTTIME(4),STARTTIME(5),STARTTIME(6))+ti/86400;

Recef = R;

for i = 1:length(ti)
    Recef(:,i) = R3(gstime(jd(i)))*R(:,i);
end

Rsun = R;

for i = 1:length(Rsun)
    Rsun(:,i) = sun(jday(STARTTIME(1),STARTTIME(2),STARTTIME(3),STARTTIME(4),STARTTIME(5),STARTTIME(6))+ti(i)/86400);
end

Rsun = Rsun*au1;

rsun = Rsun-R;
rsuntsw=rsun;
for i = 1:length(rsun)
    rsuntsw(:,i) = xyz2tsw(rsun(:,i),Om(i),in(i),u(i));
end

theta = zeros(length(ti),1);
phi = zeros(length(ti),1);

for i = 1:length(rsuntsw)
    theta(i) = acos(rsuntsw(3,i) / norm(rsuntsw(:,i)));
    phi(i) = atan2(rsuntsw(2,i),rsuntsw(1,i));
end

plot(ti,[theta,phi]/deg,'LineWidth',2);grid;

zer = zeros(length(ti),5);

MAM = [zer,ti,jd,Om,in,u,R',Recef',Rsun',theta,phi];

end

