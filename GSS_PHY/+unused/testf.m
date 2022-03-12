function [ MF ] = testf( N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

p = zeros(N,1);
l1 = zeros(N,1);
l2 = zeros(N,1);
Om = zeros(N,1);
in = zeros(N,1);
u = zeros(N,1);
m = zeros(N,1);
jd = zeros(N,1);
solar = zeros(N,3);

fgtsw = zeros(N,3);
fatmtsw = zeros(N,3);
fsunmoontsw = zeros(N,3);
femptsw = zeros(N,3);
ftsw = zeros(N,3);

for i = 1:N
    clc
    fprintf('%3.2f%% COMPLETE;\n', (i-1)/N*100);
    p(i) = 6800 + 500*rand;
    l1(i) = 1E-05 * rand;
    l2(i) = 1E-05 * rand;
    Om(i) = 2*pi*rand;
    in(i) = pi*rand;
    u(i) = 2*pi*rand;
    m(i) = 100 + 900*rand;
    jd(i) = 2457389 + 365 * rand;
    [solar(i,1), solar(i,2), solar(i,3)] = getsolar(jd(i),3);
    if i == 1
        p(i) = 600+6378.14;
        l1(i) = 0;
        l2(i) = 0;
        Om(i) = 10.9*pi/180;
        in(i) = 97.792*pi/180;
        u(i) = pi/4;
        m(i) = 100 + 900*rand;
        jd(i) = 2.457832936111111e+06;       
        solar(i,1) = 115;
        solar(i,2) = 125;
        solar(i,3) = 3;
    end
    [ fgtsw(i,:), fatmtsw(i,:), fsunmoontsw(i,:), femptsw(i,:), ftsw(i,:) ] = test1f( p(i), l1(i), l2(i), Om(i), in(i), u(i), m(i), jd(i),  solar(i,2), solar(i,1), solar(i,3) );
    
end

clc;
fprintf('TRANSFORMATION COMPLETE;\n');

MF = [p, l1, l2, Om, in, u, m, jd, solar, zeros(N,1),fgtsw, fatmtsw, fsunmoontsw, femptsw, ftsw];

end

