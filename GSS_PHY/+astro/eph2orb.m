function [ORBIT_TESSERACT] = eph2orb(eph_cube, ephcart_cube)
%EPH2ORB Summary of this function goes here
%   Detailed explanation goes here

Nsat = numel(eph_cube(:,1,1));
Nt = numel(eph_cube(1,:,1));
par = 0:360; %DEGREES

ORBIT_TESSERACT = zeros(Nsat,Nt,361,3);

if nargin < 2
    for i = 1:Nsat
        [R,V] = math.randv(eph_cube(i,:,1),eph_cube(i,:,2), ...
            eph_cube(i,:,5),eph_cube(i,:,4),eph_cube(i,:,3), ...
            eph_cube(i,:,6)-eph_cube(i,:,3));
        ephcart_cube(i,:,1:3) = R';
        ephcart_cube(i,:,4:6) = V';
    end
elseif size(eph_cube) ~= size(ephcart_cube)
    error('u mad bro?');
end

%% IMPL

pplane = eph_cube(:,:,1) .* (1 - eph_cube(:,:,2).^2); % XY plane of focal parameter
nuplane = eph_cube(:,:,6) - eph_cube(:,:,3);

for i = 1:Nsat
    for j = 1:Nt
        
    end
end

end

