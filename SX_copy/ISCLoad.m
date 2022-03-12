function [ DEPOTS, USAGE, isolated, chosenHubs ] = ISCLoad( CONN,cfg, depots, truetraffic )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pl = cfg(1); spp = cfg(2);
satcou = pl * spp;

truetraffic = truetraffic*8/1024^2;

USAGE = zeros(satcou);

if isscalar(depots) && depots < 0
    gscount = -depots;
else
    gscount = numel(depots);
end

if isempty(depots)
    DEPOTS = ceil(satcou*rand(1,gscount));
    DEPOTS = orderize(DEPOTS);
else
    DEPOTS = depots;
    DEPOTS = orderize(DEPOTS);
end

DC = length(DEPOTS);

chosenHubs = zeros(satcou);

Lmin = Inf;
isolated = [];

fprintf(['Load analysis begins >>> ',datestr(clock,0), '\n']);

for q = 1:satcou
    fprintf('-> sat #%3.0f ->',q);
    tic
    if any(DEPOTS == q)
       fprintf(' ...IS A DEPOT!\n'); 
       continue
    end
    ustemp = zeros(satcou);
    Lmin = Inf;
    if DC <= 12
        DCmod = DC;
        pldepots = DEPOTS;
    else
        pldepots  = findPlausibleHubs( q, cfg, DEPOTS, 0.5 );
        DCmod = length (pldepots);
    end
    for w = 1:DCmod
        %fprintf('...route to #%g\n',DEPOTS(w));
        try
            fprintf(' %g',pldepots(w));
            [PA,US] = findroute(q,pldepots(w),CONN,cfg);
            if (length(PA) < Lmin) || (length(PA) == Lmin && rand > 0.5)
                Lmin = length(PA);
                TRAFFIC = truetraffic(q);
                ustemp = US.*TRAFFIC;
                chosenHubs(q) = pldepots(w);
            end
        catch
            fprintf('WARNING: depot sat #%g appears to be unreachable from sat #%g!\n',pldepots(w),q);
        end %endtry
    end %endfor w
    if Lmin == Inf
        isolated = [isolated,q];
    end
    USAGE = USAGE + ustemp;
    T=toc;
    fprintf(' || %4.2f sec.\n',T);
end %endfor q


end

