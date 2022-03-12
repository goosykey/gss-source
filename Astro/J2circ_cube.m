function [eph, ephcart] = J2circ_cube(ini0,t)
%PERFORMS J2 circular perturbation for N satellites
%   
%   Initial release: 24.05.2018 by A. Kharlan
%       Returns a 3D array, thus the name xD
%   
%   INPUT:
%       ini0: Nsatx6 - keplerian a-e-om-Om-in-u intial parameters
%       t   : time points in seconds (1xNt)
%       
%   Output:
%       eph : NsatxNtx6 - keplerian for every t and every sat
%       ephcart : NsatxNtx6 - cartesian ECI for every t and every sat


Nsat = numel(ini0(:,1));

[ eph ] = astro.j2circ( ini0, t );

ephcart = eph .* 0;

for i = 1:Nsat
    [R,V] = math.randv(eph(i,:,1),eph(i,:,2),eph(i,:,5),eph(i,:,4),eph(i,:,3),eph(i,:,6)-eph(i,:,3));
    ephcart(i,:,1:3) = R';
    ephcart(i,:,4:6) = V';
end

end

