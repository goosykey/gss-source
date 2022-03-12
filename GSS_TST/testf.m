function [T, v] = testf(h)
%TESTF найти период
%   ckugckukufwkc

%% CONST

R0 = 6378.137;

mu = 398600.44;

%% IMPL

a = tst2(R0,h);

T = 2*pi*sqrt(a.^3./mu);

v = sqrt(mu/a);

end


function [a] = tst2(R0,h)

a = R0 + h;

end