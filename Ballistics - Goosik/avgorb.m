function [ yavg ] = avgorb( hrperiod, t,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if hrperiod == 0
    avgperiod = t(length(t));
else avgperiod = hrperiod*3600;
end

L = length(y);
y1 = y;
cou = 1;
num = 0;
curbeg = 1;
sum = y(1,:);
sum = sum-sum;

for i = 1:L
    if t(i) > avgperiod * cou %if exceeds
        sum = sum ./ num;
        for j = curbeg:(curbeg + num-1)
            y1(j,:) = sum;
        end
        curbeg = i;
        sum = y(i,:);
        num = 1;
        cou = cou + 1;
    else % = doesn't exceed
        sum = sum + y(i,:);
        num = num + 1;
        if t(i) == avgperiod * cou || i == L %if fin or break
            sum = sum ./num;
            for j = curbeg:(curbeg + num-1)
                y1(j,:) = sum;
            end
            curbeg = i + 1;
            sum= sum - sum;
            num = 0;
            cou = cou + 1;
        end %endif fin or break
    end %endif exceeds
    
end %i cycle

yavg = y1;

end

