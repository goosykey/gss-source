function [ LATS, AREAS ] = vis_regions_go( ini0, jd0, period, pla, spp )
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

mapq = 'low';
image_file = ['EM_' mapq '.jpg'];
cdata = imread(image_file); % Load Earth image for texture map

for i = 1:length(t)
    cla;
    PRIMARY_STATE = eph0(i,:);
    [ ~, CART510, ~ ] = sat.seedconstel( PRIMARY_STATE, jt(i), pla, spp );
    [ ~, latsat, area ] = vis_regions1( CART510, jt(i), cdata );
    add3d.majorcities(jt(i), true);
    add3d.sunlight(jt(i));
    add3d.moonpoint(jt(i));
    
    LATS(i) = latsat;
    AREAS(i) = area;
    pause(0.01);
end


end

