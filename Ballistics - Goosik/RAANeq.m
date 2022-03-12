function [ Om02, Omk2, t02, tk2 ] = RAANeq( Omega00, t01, Omega02, driftph, driftop, Dtact, deltaOmega )
%GET ACTIVE STAGE DURATION
%   Solve linear equations for Om02, Omk2, t02, tk2 

MATR = [1 0 -driftph 0
    0 1 0 -driftop
    0 0 -1 1
    -1 1 0 0];
RIGH = [Omega00 - driftph*t01; Omega02 - driftop*t01; Dtact; deltaOmega];

XVEC = MATR \ RIGH;
Om02 = XVEC(1);
Omk2 = XVEC(2);
t02 = XVEC(3);
tk2 = XVEC(4);

end

