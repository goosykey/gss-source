function [ rrelgnd, vrelgnd, vrelgndrad, visible ] = relgnd( t, rxyz, vxyz, fi, la )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

L = length (t);

vgnd = zeros(3,L);
rgnd = zeros(3,L);
visible = zeros(1,L);
vrelgndrad = zeros(1,L);

om0=7.292114e-5;
R0 = 6378.137;

global STARTTIME;

ST0 = STARTTIME;
jd0 = jday(ST0(1), ST0(2), ST0(3), ST0(4), ST0(5), ST0(6));

for i = 1:L
    jd = jd0 + t(i)/86400;
    theta = gstime(jd) + la;
    rgnd(:,i) = [R0*cos(fi)*cos(theta); R0*cos(fi)*sin(theta); R0*sin(fi)];
    vgnd(:,i) = [-om0*R0*cos(fi)*sin(theta); om0*R0*cos(fi)*cos(theta);0]; 
end

rrelgnd = rxyz-rgnd;
vrelgnd = vxyz-vgnd;

for i = 1:L
   vrelgndrad(i) = rgnd(:,i)'*vrelgnd(:,i)/norm(rgnd(:,i));
end

for i = 1:L
    if rrelgnd(:,i)'*rgnd(:,i) >=0
        visible(i) = 1;
    else
        visible(i) = nan;
    end
end

end

