function [ cdb ] = readcities2( filename, census )
%READ SHITTY CITY DATABASE
%   Detailed explanation goes here

[Country,City,~,~,Population,Latitude,Longitude] = textread(filename,'%s %s %s %s %f %f %f', -1, 'delimiter',',','headerlines',1,'emptyvalue',0);
%fid = fopen(filename);

%C = textscan(fid,'%s %s %s %s %f %f %f', -1, 'delimiter',',','headerlines',1);

%fclose(fid);

%cdb = C;
cdb = [Country, City, num2cell([Population, Latitude, Longitude])];

end

