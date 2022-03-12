function [gscart] = gs_cube(GS, t, jd0)
%Returns ground station coordinates in 3D array
%   
%   Initial release: 24.05.2018 by A. Kharlan
%       Returns a 3D array, thus the name xD
%   
%   INPUT:
%       ini0: Ngs x 2 - LAT,LON of ground stations
%       t   : time points in seconds (1xNt)
%       jd0 : julian date of t=0, scalar
%       
%   Output:
%       gscart : NgsxNtx3 - cartesian ECI for every t and every GS

R0 = 6378.137;
deg = pi/180;

Ngs = numel(GS(:,1));
jt = jd0 + t/86400;
Nt = numel(t);

gscart = zeros(Ngs,Nt,3);

for i = 1:Ngs
    X = R0 * cosd(GS(i,1)) .* cos(GS(i,2)*deg + astro.gstime(jt));
    Y = R0 * cosd(GS(i,1)) .* sin(GS(i,2)*deg + astro.gstime(jt));
    Z = R0 * sind(GS(i,1)) + jt*0;
    gscart(i,:,:) =[X; Y; Z]';
end

end

