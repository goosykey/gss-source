function [ output_args ] = plotnet( SATCOOR, CONN, cfg, plotbrokenlinks )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin<4
    plotbrokenlinks = false;
end

pl = cfg(1); spp = cfg(2);
satcou = pl * spp;
SATCOOR1 = SATCOOR;
SATCOOR1(spp:spp:satcou,2) = SATCOOR1(spp:spp:satcou,2)-360;

CONN_1 = real(CONN).*heaviside(real(CONN)) + imag(CONN).*heaviside(imag(CONN))*1i;
CONN_0 = real(CONN).*heaviside(-real(CONN)) + imag(CONN).*heaviside(-imag(CONN))*1i;

figure('Name','Plot','NumberTitle','on');
plot(SATCOOR(:,1),SATCOOR(:,2), 'o');hold on;
gplot(real(CONN_1),SATCOOR, '-g');
gplot(imag(CONN_1),SATCOOR1,'-g');
if plotbrokenlinks
    gplot(real(CONN_0),SATCOOR, ':k');
    gplot(imag(CONN_0),SATCOOR1,':k');
end
axis square;




end

