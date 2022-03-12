function [ Plmsinfi ] = legendre_my( l, m, fi )
%Universal func. for Legendre functions
%   YalinyCPB p.100

Q = [0.5 0 0 0 3 0 -1
        3 1 0 0 0 1 0
        3 2 0 0 0 0 1
        0.5 0 0 5 0 -3 0
        0.5 1 0 0 15 0 -3
        15 2 0 0 0 1 0
        15 3 0 0 0 0 1
        0.125 0 35 0 -30 0 3
        2.5 1 0 7 0 -3 0
        7.5 2 0 0 7 0 -1
        105 3 0 0 0 1 0
        105 4 0 0 0 0 1
        0 0 0 0 0 0 0];
    
if l == 2
    N = 1 + m;
elseif l == 3
    N = 4 + m;
elseif l == 4
    N = 8 + m;
end

if m > l
    Plmsinfi = 0;
else
    Plmsinfi = Q(N,1)*(cos(fi))^Q(N,2) * (Q(N,3)*(sin(fi))^4 + Q(N,4)*(sin(fi))^3 + Q(N,5)*(sin(fi))^2 + Q(N,6)*sin(fi) + Q(N,7));
end

    
end

