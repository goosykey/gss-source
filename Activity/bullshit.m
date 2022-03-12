
angles = [];
lats = [];

for poop = 0:10:90
    ini0(6) = poop*deg;
    [ STATE510, CART510, LATLON510 ] = getstate510( ini0, jd0, netdata );
    XYZ = CART510(241,1:3);  % SATELLITE CARTESIAN POSITION ECI
    Vxyz = CART510(241,4:6); % SATELLITE CARTESIAN VELOCITY ECI
    LATLON = [asin(XYZ(3)/norm(XYZ)) atan2(XYZ(2),XYZ(1))] / deg; % SUBSAT POINT
    EAST = [-sind(LATLON(2)), cosd(LATLON(2)), 0]; % EAST VECTOR @ SUBSAT
    Vsubsat = Vxyz * R0 / norm(XYZ); % SUBSAT VELOCITY
    ANGLE = acos(dot(Vsubsat, EAST) / (norm(Vsubsat) * norm(EAST))) / deg;
    angles(poop/10+1) = ANGLE;
    lats(poop/10+1) = LATLON(1);
end

angles = angles';
lats = lats';
