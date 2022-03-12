function [ weatfull, weatbool ] = getclouds( weatfull, Omask, mapdata )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[M,N] =  size (weatfull);

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

lats = latlim(1):-latstep:latlim(2);
lons = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(lats,lons);

Omask = boolean(imresize(double(Omask),[M N],'bilinear'));

latmap = imresize(latmap,[M N],'bilinear');
lonmap = imresize(lonmap,[M N],'bilinear');

%% NOMINAL WEATHER EFFECTS

weat1 = rand(M,N)-0.5; % all random -0.5..0.5

% LATITUDE EFFECTS

aa = latmap(:,1)';
weat2_0 = 0.4*exp(-((aa)/10).^2) + 0.025*randn(1,numel(aa)); % cloud belt - rainforest
weat2_1 = 0.4*exp(-((aa)/30).^2) + 0.025*randn(1,numel(aa)); % cloud belt - ocean

weat2_2 = -0.3*exp(-((abs(aa)-23.19)/20).^2) + 0.025*randn(1,numel(aa)); % drought - tropical

weat2_0 = ~Omask .*repmat(weat2_0',[1,N]);
weat2_1 = Omask  .*repmat(weat2_1',[1,N]);
weat2_2 = ~Omask .*repmat(weat2_2',[1,N]);

weat2 = weat2_0 + weat2_1 + weat2_2;

% TERAIN TYPE EFFECTS

weat3 = ~Omask .* rand(M,N).*0.1 .* (-1); % slightly improve inland weather
weat4 =  Omask .* rand(M,N).*0.1; % slightly impair ocean weather

weatnom = weat1 + weat2 + weat3 + weat4;
weatnom(weatnom < -1) = -1;
weatnom(weatnom >  1) =  1;

%% WEATHER DEVIATIONS

delta = weatfull-weatnom;

weatfull = weatfull - 0.6* delta + 0.25*randn(M,N);
weatfull(weatfull < -1) = -1;
weatfull(weatfull >  1) =  1;
%weatfull = weatnom;


weatbool = boolean(heaviside(weatfull-0.1));

end

