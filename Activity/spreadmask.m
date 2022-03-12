function [ maskNEW ] = spreadmask( mask, R, C )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

maskNEW = int8(zeros(R,C));
[r,c] = size(mask);

for i = 1:R
    for j = 1:C
        maskNEW(i,j) = mask(ceil(i*r/R),ceil(j*c/C));
    end
end

end

