function [ ANS ] = timetime( S, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

yr1=S(1);
mon1=S(2);
day1=S(3);
hr1=S(4);
min1=S(5);
sec1=S(6)+t;

while sec1 >= 60
    sec1 = sec1-60;
    min1 = min1 + 1;
end

while min1 >= 60
    min1 = min1-60;
    hr1 = hr1 + 1;
end

while hr1 >= 24
   hr1 = hr1 - 24;
   day1 = day1 +1;
   if ismember(mon1,[1,3,5,7,8,10,12]) && day1 == 32
       day1 = 1;
       mon1 = mon1+1;
   elseif ismember(mon1,[4,6,9,11]) && day1 == 31
       day1 = 1;
       mon1 = mon1+1;
   elseif mon1==2 &&day1==(30-sign(mod(yr1,4)))
        day1 = 1;
       mon1 = mon1+1;
   end
end

while mon1 > 12
    mon1 = mon1-12;
    yr1 = yr1 + 1;
end

ANS = [yr1, mon1, day1, hr1, min1, sec1];

end

