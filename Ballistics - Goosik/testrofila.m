function [ MRFL ] = testrofila( N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

rofila = zeros(N,3);
jd = zeros(N,1);
dUro = zeros(N,1);
dUfi = zeros(N,1);
dUla = zeros(N,1);

a = 6500;
b = 7500;

for i = 1:N
    rofila(i,1) = a + (b-a)*rand;
    rofila(i,2) = -pi/2 + pi*rand;
    rofila(i,3) = 2*pi*rand;
    jd(i) = 2451545.0 + 12000 * rand;
    [ dUro(i), dUfi(i), dUla(i) ] = rofila1test( rofila(i,1), rofila(i,2), rofila(i,3));
end

MRFL = [rofila,jd,dUro,dUfi,dUla];
fprintf('TRANSFORMATION COMPLETE;\n');

end

