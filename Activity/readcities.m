function [ cdb ] = readcities( filename, census )
%READ GEONAME.ORG CITY DATABASE
%   INPUT:
%   filename : string file name, e.g. city1000.txt
%   census   : removes all entries with population less then this number
%   
%   OUTPUT:
%   cdb      : cell array city database

fprintf('Reading file data...\n');
% [geonameid, name, asciiname, alternatenames, latitude, longitude, featclass,...
%     featcode, countrycode, cc2, a1c, a2c, a3c, a4c, population, elevation,...
%     dem, timezone, modifdate] = textread(filename,...
%     '%f %s %s %s %f %f %s %s %s %s %s %s %s %s %f %f %f %s %s',...
%     -1, 'delimiter','\t','headerlines',0);

[~, ~, asciiname, ~, latitude, longitude, ~,...
    ~, countrycode, ~, ~, ~, ~, ~, population, ~,...
    ~, ~, ~] = textread(filename,...
    '%f %s %s %s %f %f %s %s %s %s %s %s %s %s %f %f %f %s %s',...
    -1, 'delimiter','\t','headerlines',0);


%fid = fopen(filename);

%C = textscan(fid,'%s %s %s %s %f %f %f', -1, 'delimiter',',','headerlines',1);

%fclose(fid);

%cdb = C;
%cdb = {Country, [Latitude, Longitude]};

%% COMPUTE ZONAL MAP

fprintf('Computing zonal map...\n');
[ Rurb, Rrur ] = twozones( population );

Rurb(Rurb(:)<0) = 0;

cdb = [countrycode, asciiname, num2cell([population, latitude, longitude, Rurb, Rrur])];

%% APPLY CENSUS IF SPECIFIED

if nargin > 1
    fprintf('Applying census...\n');
    L1 = length(cdb(:,3));
    cdb(cell2mat(cdb(:,3))<census,:) = [];
    L2 = length(cdb(:,3));
    fprintf('%g out of %g entries left (%g removed).\n',L2,L1,(L1-L2));
end

fprintf('Done.\n');
end

