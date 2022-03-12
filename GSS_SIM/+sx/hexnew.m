function [ hexacx, hexacy, hexasx, hexasy ] = hexnew( lonc, latc, Rcov, a )
%HEXNEW Summary of this function goes here
%   Detailed explanation goes here

b = a*sqrt(3)/2;
c = a/2;

x0 = 0;
y0 = -Rcov;

hexacx = []; hexacy = [];
hexasx = []; hexasy = [];
cou = 0;

while y0 < Rcov
    [xx, yy] = genhex (x0, y0, a, b, c);
    hexacx = [hexacx;x0]; hexacy = [hexacy;y0];
    hexasx = [hexasx;xx]; hexasy = [hexasy;yy];
    
    [xx, yy] = genhex (x0+a+c, y0+b, a, b, c);
    hexacx = [hexacx;x0+a+c]; hexacy = [hexacy;y0+b];
    hexasx = [hexasx;xx]; hexasy = [hexasy;yy];
    
    [xx, yy] = genhex (x0-a-c, y0+b, a, b, c);
    hexacx = [hexacx;x0-a-c]; hexacy = [hexacy;y0+b];
    hexasx = [hexasx;xx]; hexasy = [hexasy;yy];
    
    cou = 1;
    
    while sqrt((x0+3*a*cou)^2 + y0^2) < Rcov
        [xx, yy] = genhex (x0+3*a*cou, y0, a, b, c); % + axial
        hexacx = [hexacx;x0+3*a*cou]; hexacy = [hexacy;y0];
        hexasx = [hexasx;xx]; hexasy = [hexasy;yy];
        
        [xx, yy] = genhex (x0+3*a*cou+a+c, y0+b, a, b, c); % + up
        hexacx = [hexacx;x0+3*a*cou+a+c]; hexacy = [hexacy;y0+b];
        hexasx = [hexasx;xx]; hexasy = [hexasy;yy];
        
        [xx, yy] = genhex (x0-3*a*cou, y0, a, b, c); % - axial
        hexacx = [hexacx;x0-3*a*cou]; hexacy = [hexacy;y0];
        hexasx = [hexasx;xx]; hexasy = [hexasy;yy];
        
        [xx, yy] = genhex (x0-3*a*cou-a-c, y0+b, a, b, c); % - up
        hexacx = [hexacx;x0-3*a*cou-a-c]; hexacy = [hexacy;y0+b];
        hexasx = [hexasx;xx]; hexasy = [hexasy;yy];
        
        cou = cou+1;
    end
    
    y0 = y0 + 2*b;
    
end

centercoords = [hexacx,hexacy];
radii = sqrt(sum(centercoords.^2,2));
excess = radii > Rcov;

hexacx(excess) = [];
hexacy(excess) = [];
hexasx(excess,:) = [];
hexasy(excess,:) = [];

deg1 = 6378.14 * 2*pi / 360;
scale = 1/deg1;

hexacy = hexacy * scale + latc;
hexasy = hexasy * scale + latc;

hexacx = hexacx * scale ./ cosd(hexacy) + lonc;
hexasx = hexasx * scale ./ cosd(hexasy) + lonc;

end

function [xx, yy] = genhex (x, y, a, b, c)

xx = [-a -c c a c -c];
yy = [0 b b 0 -b -b];

xx = xx+x;
yy = yy+y;

end


