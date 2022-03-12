function all_coverage( CARTN, alphacrit, circleformat, conecolor )
%VIS_COVERAGE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    conecolor = 'c';
    if nargin < 3
        circleformat = ':r';
    end
end

add3d.plotcov(CARTN, alphacrit, circleformat);
add3d.plotcone(CARTN,alphacrit, conecolor);


end

