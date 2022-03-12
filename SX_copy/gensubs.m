function [ map1, map2, map3 ] = gensubs( map, mapdata, MATRIX, medianstack, sigma)
%GENERATE OCEAN DWELLERS MASK
%   Detailed explanation goes here

% MATRIX: 
% 
%   TYPE1P TYPE1I TYPE1S
%   TYPE2P TYPE2I TYPE2S
%   TYPE3P TYPE3I TYPE3S

[M,N] = size(map);
MAP3D = zeros(M,N,3);

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

if latlim(1)>latlim(2)
    latstep = -latstep;
end

LATS = latlim(1):latstep:latlim(2);
LONS = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(LATS,LONS);

%% FORM OCEAN MASK

tic;
fprintf('Forming ocean mask...\n');
Omask = oceanALL(latmap,lonmap,10);
T = toc;
fprintf('Mask formed in %3.2f sec.\n',T);

Landmask = not(Omask);

for q = 1:3
    
    peoplesubs = MATRIX(q,1);
    landsubs = MATRIX(q,2);
    straysubs = MATRIX(q,3);
    
    
    %% GENERATE PEOPLE-RELATED SUBS
    
    maincoef = peoplesubs/sum(sum(map));
    map_peoplesubs = map .* maincoef;
    
%     %% GENERATE INVERSE PEOPLE-RELATED SUBS
%     
%     maxmax = max(max(map));
%     map_reduced = map./maxmax;
%     
%     inv_map_reduced = (ones(M,N) - map_reduced).*Landmask;
%     maincoef = invsubs/sum(sum(inv_map_reduced));
%     map_invsubs = inv_map_reduced .* maincoef;

    %% GENERATE LAND RANDOM STRAY SUBS
    
    stacksnum = round(landsubs/medianstack);
    
    fprintf('%3.0f land stacks to be generated... \n',stacksnum);
    
    latrand = ceil(rand(stacksnum,1)*M);
    lonrand = ceil(rand(stacksnum,1)*N);
    coorand = [latrand,lonrand];
    
    map_landsubs = zeros(M,N);
    
    for i = 1:stacksnum
        stack = sigma*randn*medianstack + medianstack;
        if stack<0
            continue
        end
        map_landsubs(coorand(i,1),coorand(i,2)) = stack;
    end
    
    map_landsubs = map_landsubs.*Landmask;
    maincoef = landsubs / sum(sum(map_landsubs));
    
    map_landsubs = map_landsubs.*maincoef;
    
    
    %% GENERATE UNINFORMLY RANDOM STRAY SUBS
    
    stacksnum = round(straysubs/medianstack);
    
    fprintf('%3.0f stray stacks to be generated... \n',stacksnum);
    
    latrand = ceil(rand(stacksnum,1)*M);
    lonrand = ceil(rand(stacksnum,1)*N);
    coorand = [latrand,lonrand];
    
    map_straysubs = zeros(M,N);
    
    for i = 1:stacksnum
        stack = sigma*randn*medianstack + medianstack;
        if stack<0
            continue
        end
        map_straysubs(coorand(i,1),coorand(i,2)) = stack;
    end
    
    %% Finally
    
    %MAP3D(:,:,q) = map_peoplesubs + map_landsubs + map_straysubs;
    MAP3D(:,:,q) = map_straysubs;
    
    
    
end

map1 = MAP3D(:,:,1);
map2 = MAP3D(:,:,2);
map3 = MAP3D(:,:,3);


end

