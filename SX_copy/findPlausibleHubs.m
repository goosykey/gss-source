function [ pldepots ] = findPlausibleHubs( sat, cfg, DEPOTS, FACTOR )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

L = length(DEPOTS);

depotable = ([DEPOTS',zeros(L,1)]);

for q = 1:L
    depotable(q,2) = satdistance(sat,depotable(q,1),cfg);
end

while length(depotable(:,1)) >= L*FACTOR
    [~,maxdi] = max(depotable(:,2));
    depotable(maxdi,:) = [];
end

pldepots = depotable(:,1)';

end

