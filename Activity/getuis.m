function [ VUI, DUI, SUI, PEN, SPEN, CountryNo ] = getuis( isocode, uidb )
%GET Usage Intensity from DATABASE
%   Detailed explanation goes here

codes = uidb(:,1);
uis = cell2mat(uidb(:,2:6));
VUI = 0;
DUI = 0;
SUI = 0;
PEN = 0.005;
SPEN = 2e-4;
CountryNo = 9;

if strcmp(isocode, 'XK') % KOSOVO IS NOT ON ISO DATABASE
    isocode = 'RS'; % KOSOVO IS SERBIA!
end

if strcmp(isocode, 'WW')
    VUI = 50;
    DUI = 50;
    SUI = 50;
    return
end

for i = 1:length(codes)
    if codes{i} == isocode
        VUI = uis(i,1);
        DUI = uis(i,2);
        SUI = uis(i,3);
        PEN = uis(i,4);
        SPEN= uis(i,5);
        CountryNo = i;
        break
    end
end


end

