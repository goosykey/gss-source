function [ mini ] = minopenfscore( openSet, fScore )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L = length(fScore);


for q = 1:L
    if ~any(openSet==q)
        fScore(q)=Inf;
    end
end

[~,mini] = min(fScore);

end

