function [ neig ] = neighbours( CONN, satno )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

neig = [];

for i = 1:length(CONN(1,:))
    if CONN(satno,i) ~= 0
        neig = [neig,i];
    end
end

end

