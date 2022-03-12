function [ totalpath ] = getpath( cameFrom, current )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

totalpath = current;
%while any (cameFrom==current)
while current ~= 0
    current = cameFrom(current);
    totalpath = [totalpath,current];
end

L = length(totalpath);

totalpath = totalpath(L-1:-1:1);


end

