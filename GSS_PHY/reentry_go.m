iniSK
PREC = 0;
[ tact,yact,tfall,yfall ] = reentry( jd0, Y );
y = [yact]; t = [tact];
p = y(:,1);
l1 = y(:,2);
l2 = y(:,3);
e = sqrt(l1.^2 + l2.^2);
a = p./(1-e.^2);
om = atan2(l2,l1);
Om = y(:,4);
in = y(:,5);
u = y(:,6);
[R,~] = math.randv(a,e,in,Om,om,u-om);

figure
plot(t/86400-jd0,sqrt(sum(R.^2)));
xlabel('time, days'); ylabel('radius, km');
figure
plot(t/86400-jd0,e);
xlabel('time, days'); ylabel('eccentricity');
m = y(:,7);
figure
plot(t/86400-jd0,m);
xlabel('time, days'); ylabel('mnass, kg');