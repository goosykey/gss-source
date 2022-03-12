function [ peoplenew ] = convertmesh( map, mapdata, F, V )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.14;
epsdeg = 0.5;

peoplenew = F(:,1)*0;

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

if latlim(1)>latlim(2)
    latstep = -latstep;
end

LATS = latlim(1):latstep:latlim(2);
LONS = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(LATS,LONS);

latL = length(LATS);
lonL = length(LONS);

latmapnew = latmap + latstep/2;
lonmapnew = lonmap + lonstep/2;
mapnew = reshape(map,[numel(map),1]);

pointtable = [reshape(latmapnew,[numel(latmapnew),1]),...
    reshape(lonmapnew,[numel(lonmapnew),1]),...
    (1:numel(latmapnew))'];

Vgeo = [asind(V(:,3)) , atan2d(V(:,2),V(:,1))];

Nfaces = length(F(:,1));

fprintf('Go cycle\n');

parpool(3);
fprintf('Parallel latitude scanning begins... \n');

parfor i = 1:Nfaces
    
%     if mod(i,5) ~= 0
%         continue
%     end
    tic
    xv = Vgeo(F(i,(~isnan(F(i,:)))),2);
    yv = Vgeo(F(i,(~isnan(F(i,:)))),1);
    
    yc = sum(yv)/numel(yv); xc = sum(xv)/numel(xv);
    pthere = trunctrunc(yc, xc, epsdeg, pointtable);
    toto1 = toc;   
    
    tic
    inporaw = inpolygon(pthere(:,2),pthere(:,1),xv,yv);
    inpo = pthere(inporaw,3);
    peoplenew(i) = sum(mapnew(inpo));

    toto = toc;
    if mod(i,10) == 0
        fprintf('step %g ; trunc in %3.2f sec ; total %3.2f sec ; %g of %g lie ; %3.2f%% done\n',...
            i,toto,toto+toto1,sum(inporaw(inporaw~=0)),length(pthere(:,1)),i/Nfaces*100);
    end
end

delete(gcp('nocreate'));
fprintf('\n');

end

function [ptnew] = trunctrunc (latappr, lonappr, eps, ptable)

tic

% ptable(abs(ptable(:,1)-latappr) > eps,:) = [];
% ptable(abs(ptable(:,2)-lonappr) > eps*cosd(latappr),:) = [];

ptable = ptable(abs(ptable(:,1)-latappr) < eps,:);
ptable = ptable(abs(ptable(:,2)-lonappr) < eps*cosd(latappr),:);

ptnew = ptable;

toto1 = toc;
    
 
    
end








