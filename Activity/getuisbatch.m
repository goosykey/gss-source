function [ VUI, DUI, SUI, PEN, SPEN, CountryNo ] = getuisbatch( isocode, uidb )
%GET Usage Intensity from DATABASE
%   ISOCODE IS A Nx1 CELL ARRAY

codes = uidb(:,1);
uis = cell2mat(uidb(:,2:6));

isocode(ismember(isocode,'XK')) = {'RS'}; % KOSOVO IS SERBIA!

[~,ism2] = ismember(isocode,codes);
CountryNo = ism2;
VUI = uis(ism2,1);
DUI = uis(ism2,2);
SUI = uis(ism2,3);
PEN = uis(ism2,4);
SPEN= uis(ism2,5);


end

