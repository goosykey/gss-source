function [ output_args ] = replot( rcv, mapdata )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



rcv = single(rcv);

rcv (~rcv) = nan;

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

lats = latlim(1):-latstep:latlim(2);
lons = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(lats,lons);

for i = 1:numel(rcv)
    if isnan(rcv(i))
        lat = latmap(i);
        lon = lonmap(i);
        rectangle('Position',[lat lon latstep lonstep],'FaceColor',[0 1 0.5]);
    end
end

end

