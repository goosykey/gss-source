function epheplot2d( jd0,t,y )
%EPHEPLOT2D Summary of this function goes here
%   Detailed explanation goes here

mu = 398600.44;

p = y(:,1);
l1 = y(:,2);
l2 = y(:,3);
e = sqrt(l1.^2 + l2.^2);
a = p./(1-e.^2);
om = atan2(l2,l1);
Om = y(:,4);
in = y(:,5);
u = y(:,6);

Pkepl = 2*pi * a(1) * sqrt(a(1) / mu); % keplerian period

jt = jd0 + t/86400;

[R, ~] = math.randv(a,e,in,Om,om,u-om);

A = math.R3v(astro.gstime(jt));

n = numel(jt);
Recf = zeros(3,n);

tete = cell(0);

[~, day, ~] = astro.gdate (jt);
hour = (day - floor(day))*24;

for i = 1:n
    Recf(:,i) = A(:,:,i) * R(:,i);
    if mod(t(i),1800) == 0
        tete = [tete, sprintf('%3.1f',hour(i))];
    end
end

Recf = Recf';

lateph = asind(Recf(:,3)./sqrt(sum(Recf.^2,2)));
loneph = atan2d(Recf(:,2),Recf(:,1));
deltalon = loneph(2:end) - loneph(1:end-1);
loneph(deltalon > 200) = nan;
lateph(deltalon > 200) = nan;

figure;

plot_google_map('maptype','hybrid');

loneph1 = loneph(1:round(end/7));
loneph2 = loneph(round(end/7)+1:2*round(end/7));
loneph3 = loneph(2*round(end/7)+1:3*round(end/7));
loneph4 = loneph(3*round(end/7)+1:4*round(end/7));
loneph5 = loneph(4*round(end/7)+1:5*round(end/7));
loneph6 = loneph(5*round(end/7)+1:6*round(end/7));
loneph7 = loneph(6*round(end/7)+1:end);

lateph1 = lateph(1:round(end/7));
lateph2 = lateph(round(end/7)+1:2*round(end/7));
lateph3 = lateph(2*round(end/7)+1:3*round(end/7));
lateph4 = lateph(3*round(end/7)+1:4*round(end/7));
lateph5 = lateph(4*round(end/7)+1:5*round(end/7));
lateph6 = lateph(5*round(end/7)+1:6*round(end/7));
lateph7 = lateph(6*round(end/7)+1:end);

plot(loneph1,lateph1,'color',[1,0,0]);
plot(loneph2,lateph2,'color',[1,0.5,0]);
plot(loneph3,lateph3,'color',[1,1,0]);
plot(loneph4,lateph4,'color',[0,1,0]);
plot(loneph5,lateph5,'color',[0,1,1]);
plot(loneph6,lateph6,'color',[0,0,1]);
plot(loneph7,lateph7,'color',[1,0,1]);



plot(loneph(mod(t,1800)==0),lateph(mod(t,1800)==0),'or',...
    'markersize',6,'markerfacecolor','red');
text(loneph(mod(t,1800)==0),lateph(mod(t,1800)==0),tete,...
    'verticalalignment','bottom','color','red','fontsize',8);


end

