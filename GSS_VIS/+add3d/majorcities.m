function majorcities (jd, minor)
%MAJORCITIES - add to earth GSS vis 3d plot
%   Detailed explanation goes here

if nargin < 2
    minor = false;
end

[lat,lon,names] = textread('majorcities.txt','%f %f %s', 'delimiter',', ');

triplet_major = [0 201 87];
triplet_major = triplet_major / norm(triplet_major);

triplet_minor = [56	142	142];
triplet_minor = triplet_minor / norm(triplet_minor);

for i = 1:length(names)
    names{i}(names{i}=='_') = ' ';
end

add3d.elem.city3d(jd, lat,lon,names, 'color', triplet_major, 'markersize',4,'fontsize',7);

if minor
    [lat,lon,names] = textread('minorcities.txt','%f %f %s', 'delimiter',', ');
    
    for i = 1:length(names)
        names{i}(names{i}=='_') = ' ';
    end
    
    add3d.elem.city3d(jd, lat,lon,names, 'color', triplet_minor, 'markersize',3,'fontsize',5);
end

end

