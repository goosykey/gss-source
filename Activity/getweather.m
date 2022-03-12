function [ weatfull, weatbool ] = getweather( weatfull, Omask )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[M1,N1] =  size (weatfull);

N = 48;
M = round(M1/N1*N);


lats = 90:-1:-90;
lons = -180:179;

[latmap, ~] = ndgrid(lats,lons);
latmap = imresize(latmap,[M N],'bilinear');

Omask = boolean(heaviside(imresize(double(Omask),[M N],'bilinear')-0.5));



%% NOMINAL WEATHER EFFECTS

weat1 = 0.5./cosd(latmap).*randn(M,N); % all randn

% LATITUDE EFFECTS

aa = latmap(:,1)';
weat2_0 = 0.5*exp(-((aa)/5).^2) + 0.025*randn(1,numel(aa)); % cloud belt - rainforest
weat2_1 = 0.6*exp(-((aa)/5).^2) + 0.025*randn(1,numel(aa)); % cloud belt - ocean

weat2_2 = -0.3*exp(-((abs(aa)-23.19)/120).^2) + 0.025*randn(1,numel(aa)); % drought - tropical

weat2_0 = ~Omask .*repmat(weat2_0',[1,N]);
weat2_1 = Omask  .*repmat(weat2_1',[1,N]);
weat2_2 = ~Omask .*repmat(weat2_2',[1,N]);

weat2 = weat2_0 + weat2_1 + weat2_2;

% TERAIN TYPE EFFECTS

weat3 = ~Omask .* rand(M,N).*0.05 .* (-1); % slightly improve inland weather
weat4 =  rand(M,N).*0.2; % slightly impair ALL weather

weatnom = weat1 + weat2 + weat3 + weat4;
weatnom(weatnom < -1) = -1;
weatnom(weatnom >  1) =  1;

%% WEATHER DEVIATIONS

weatfull = imresize(weatfull,[M N],'bilinear');
delta = weatfull - weatnom;
if mean(mean(abs(delta))) > 0.6
    %fprintf('!');
    weatfull = (weatnom + weatfull) /2 ;
else
    weatfull = weatfull + 0.25*randn(M,N);
end

weatfull(weatfull < -1) = -1;
weatfull(weatfull >  1) =  1;
%weatfull = weatnom;

weatfull = imresize(weatfull,[M1 N1],'bilinear');

weatbool = boolean(heaviside(weatfull-0.3));

end

