function [ output_args ] = singleframe( t, R)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

H = earthdraw(t);
L = length(R(1,:));
figure(H);

STR = ['or'; 'og'; 'ob'];

for i = 1:L
    plot3 (R(1,i)*1000, R(2,i)*1000, R(3,i)*1000, STR(i,:));
end

end

