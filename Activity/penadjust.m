function [ Pout, SPout, temp, appNsub, appSsub ] = penadjust( Pmask, SPmask, Cnummask, map, uidb )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fprintf('Adjusting penetration mask... \n');
L = length(uidb(:,1));

Pout = Pmask;
SPout = SPmask;

appNsub = map .* Pmask;
appSsub = map .* SPmask;

for i = 1:L
    fprintf ('%2.1f%% >> ',i/L*100);
    %isoc = uidb{i,1};
    meanP = uidb{i,5};
    meanSP = uidb{i,6};
    appsumP = sum(appNsub(Cnummask(:)==i));
    appsumSP = sum(appSsub(Cnummask(:)==i));
    summap = sum(map(Cnummask(:)==i));
    appmeanP = appsumP/summap;
    appmeanSP = appsumSP/summap;
    Padj = meanP/appmeanP;
    SPadj = meanSP/appmeanSP;
    Pout(Cnummask(:)==i) = Pmask(Cnummask(:)==i).*Padj;
    SPout(Cnummask(:)==i) = SPmask(Cnummask(:)==i).*SPadj;
    if i == 237
        temp = [meanP, meanSP, appsumP, appsumSP, summap];
    end
end

fprintf('Done. \n');


end

