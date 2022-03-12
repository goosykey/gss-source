function epheplot( jd0,t,y )
%EPHEPLOT Summary of this function goes here
%   Detailed explanation goes here

[ ~ ] = vis_earthdraw(jd0+t(end), 'axes', 'eci', 'prime', 'none','quality','hd');

jt= jd0+t;

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
R = R'*1e3;

plot3(R(:,1),R(:,2),R(:,3),'-c');
plot3(R(end,1),R(end,2),R(end,3),'oc', 'markerfacecolor',[0.6 0.6 0.6]);

end

