function [ reci, jd, rteme ] = testeci2teme( N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

reci = zeros(N,3);
jd = zeros(N,1);
rteme = zeros(N,3);

a = 3000;
b = 6000;

for i = 1:N
    %clc
    %fprintf('%3.2f%% COMPLETE;\n', (i-1)/N*100);
    reci(i,:) = a + (b-a)*rand(1,3);
    jd(i) = 2451545.0 + 12000 * rand;
    [reci(i,:),jd(i)]
    ttt = (jd(i) - 2451545.0  )/ 36525.0; %Julian Cent
    rteme(i,:) = eci2prec(reci(i,:)',[0;0;0],[0;0;0],ttt,106,0,'a');
end

%clc;
fprintf('TRANSFORMATION COMPLETE;\n');

end

