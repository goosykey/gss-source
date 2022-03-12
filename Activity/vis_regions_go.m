function [ LATS, AREAS ] = vis_regions_go( ini0, jd0, period, netdata )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

eph0 = J2pert(ini0,period,10);
t = eph0(:,1);
eph0 = eph0(:,2:7);
jt = jd0 + t/86400;

vis_earthdraw(jd0);

LATS = zeros(length(t),1);
AREAS = zeros(length(t),1);

global temp temp1 temp2;
temp = []; temp1 = []; temp2 = [];

for i = 1:length(t)
    cla;
    PRIMARY_STATE = eph0(i,:);
    [ ~, CART510, ~ ] = getstate510( PRIMARY_STATE, jt(i), netdata );
    [ ~, latsat, area ] = vis_regions1( CART510, jt(i) );
    LATS(i) = latsat;
    AREAS(i) = area;
    pause(0.01);
end


end

