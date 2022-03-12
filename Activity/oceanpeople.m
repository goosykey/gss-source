function [ map_oc ] = oceanpeople( mapdata, straypeople, shorepeople, medianstack, sigma, RANGEmask, OCEANmask )
%GENERATE OCEAN DWELLERS MASK
%   Detailed explanation goes here

%% INITIALIZE MAP PARAMETERS

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

lats = latlim(1):-latstep:latlim(2);
lons = lonlim(1):lonstep:lonlim(2);

%[latmap, lonmap] = ndgrid(lats,lons);

%% FORM BOOLEAN OCEAN MASK

tic;
boob = OCEANmask;

[M,N] = size(boob);

map_oc = zeros(M,N);
map_sh1 = zeros(M,N);
map_sh2 = zeros(M,N);
map_sh3 = zeros(M,N);

cell_sh = {map_sh1,map_sh2,map_sh3};

%% GENERATE UNINFORMLY RANDOM STRAY OCEAN PEOPLE

oceantoall = sum(sum(boob))/(M*N);

stacksnum = round(straypeople/medianstack/oceantoall);

fprintf('%3.0f stray stacks to be generated... \n',stacksnum*oceantoall);

latrand = ceil(rand(stacksnum,1)*M);
lonrand = ceil(rand(stacksnum,1)*N);
coorand = [latrand,lonrand];

for i = 1:stacksnum
    stack = sigma*randn*medianstack + medianstack;
    if stack<0
        continue
    end
    map_oc(coorand(i,1),coorand(i,2)) = stack;
end

map_oc = map_oc.*boob;

%% GENERATE SHORE PEOPLE

RANGEmask(RANGEmask(:)<1e-3) = 1;

stacksnum = round(shorepeople/medianstack);
fprintf('%3.0f shore stacks to be generated... \n',stacksnum);
if stacksnum >= 2000
    par = true;
    fprintf('\t too many stacks, initiating parallel calculation... \n');
else
    par = false;
end

RANGEmask = imresize(RANGEmask,[M N],'bilinear');
boobmask = boob ./ RANGEmask;

if par
    parpool(3);
    parfor I = 1:3
        for i = 1:round(stacksnum/3)
            %fprintf('\tStack # %3.0f \n',i);
            stack = sigma*randn*medianstack + medianstack;
            if stack<0
                continue
            end
            [lonrd,latrd] = pinky(lons,lats,boobmask,3);
            [~,q] = min(abs(lats-latrd));
            [~,w] = min(abs(lons-lonrd));
            cell_sh{I}(q,w) = cell_sh{I}(q,w)+stack;
        end
    end
    
    map_sh = cell_sh{1}+cell_sh{2}+cell_sh{3};
    map_oc = map_oc + map_sh;
    delete(gcp('nocreate'));
else
    for I = 1:4
        for i = 1:round(stacksnum/4)
            %fprintf('\tStack # %3.0f \n',i);
            stack = sigma*randn*medianstack + medianstack;
            if stack<0
                continue
            end
            [lonrd,latrd] = pinky(lons,lats,boobmask,3);
            [~,q] = min(abs(lats-latrd));
            [~,w] = min(abs(lons-lonrd));
            cell_sh{I}(q,w) = cell_sh{I}(q,w)+stack;
        end
    end
    map_sh = cell_sh{1}+cell_sh{2}+cell_sh{3}+cell_sh{4};
    map_sh = gather(map_sh);
    map_oc = map_oc + map_sh;
end %if par

toc;

end

