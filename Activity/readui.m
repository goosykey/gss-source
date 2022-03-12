function [ uidb ] = readui( filename )
%READ .CSV USAGE INTENSITY DATABASE
%   Detailed explanation goes here

[isocode, VUI, DUI, SUI, PEN, SPEN] = textread(filename,'%s %f %f %f %f %f', -1, 'delimiter',';','emptyvalue',10);

uidb = [isocode, num2cell([VUI, DUI, SUI, PEN, SPEN])];

end

