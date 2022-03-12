function [ val, val_ortho ] = satdistance( s1, s2, cfg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pl = cfg(1); spp = cfg(2);
satcou = pl * spp;

pno = [idivide(s1-1,int32(spp))+1 , idivide(s2-1,int32(spp))+1]; % plane num
sip = [mod(s1-1,spp)+1 , mod(s2-1,spp)+1]; % num of sat in the plane

dist_one = sqrt(double( (pno(2)-pno(1))^2 + (sip(2)-sip(1))^2 ));
dist_two = sqrt(double( (pno(2)-pno(1))^2 + (sip(2)+spp-sip(1))^2 ));
dist_three = sqrt(double( (pno(2)-pno(1))^2 + (sip(2)-spp-sip(1))^2 ));

val= min([dist_one, dist_two, dist_three]);

if nargout > 1
    ortho_one = abs(pno(2)-pno(1)) + abs(sip(2)-sip(1));
    ortho_two = abs(pno(2)-pno(1)) + abs(sip(2)+spp-sip(1));
    ortho_three = abs(pno(2)-pno(1)) + abs(sip(2)-spp-sip(1));
    val_ortho= min([ortho_one, ortho_two, ortho_three]);
end

end

