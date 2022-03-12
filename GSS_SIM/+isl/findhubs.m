function [ hubnums ] = findhubs( latlons, jd, satarray, elevcritGS )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.137;


[M,N] = size(latlons);
if N ~=3
    error('gtfo latlons');
end

[~,N1] = size(satarray);
if N1 ~=3
    error('gtfo sats');
end

hubs = [];

for i = 1:M
    z0ecf = R0*sind(latlons(i,1));
    x0ecf = R0*cosd(latlons(i,1))*cosd(latlons(i,2));
    y0ecf = R0*cosd(latlons(i,1))*sind(latlons(i,2));
    R0eci = R3(-gstime(jd))*[x0ecf y0ecf z0ecf]';
    
    [ visible, deltaangle ] = ifvisible ( R0eci, satarray, elevcritGS );
    
    deltaangle(not(visible)) = -inf;
    
    for j = 1:latlons(i,3)
        [maxda, maxi] = max(deltaangle);
        
        while any(hubs == maxi)
            if maxda < 0
                break
            end
            deltaangle(maxi) = -inf;
            [maxda, maxi] = max(deltaangle);
        end
        
        if maxda > 0
            hubs = [hubs,maxi];                
        end
        
    end
    
    
end

hubnums = orderize(hubs);

end

