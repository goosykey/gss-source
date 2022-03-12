function [ PATH, USAGE ] = findroute ( s1, s2, CONN, cfg ) 
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

pl = cfg(1); spp = cfg(2);
satcou = pl * spp;
% SATCOOR1 = SATCOOR;
% SATCOOR1(spp:spp:satcou,2) = SATCOOR1(spp:spp:satcou,2)-360;



CONN_1 = real(CONN).*heaviside(real(CONN)) + imag(CONN).*heaviside(imag(CONN))*1i;
CONN_0 = real(CONN).*heaviside(-real(CONN)) + imag(CONN).*heaviside(-imag(CONN))*1i;

USAGE = zeros(satcou);

CONN_1 = mirrorize(CONN_1);


closedSet = [];
openSet = s1;

cameFrom = zeros(satcou,1);

gScore = zeros(satcou,1); gScore(:)=Inf;
gScore(s1) = 0;

fScore = zeros(satcou,1); fScore(:)=Inf;
fScore(s1) = satdistance(s1,s2,cfg);

while ~isempty(openSet)
    %openSet
    %closedSet
    mini = minopenfscore( openSet, fScore );
    current = mini;
    if current == s2
        PATH = getpath( cameFrom, current );
    end
    openSet = openSet(openSet~=current);
    closedSet = [closedSet,current];
    closedSet = orderize(closedSet);
    neigcur = neighbours(CONN_1,current);
    for q = 1:length(neigcur)
        if any(closedSet==neigcur(q))
            continue
        end
        tent_gS = gScore(current) + 1;
        if ~any(openSet==neigcur(q))
            openSet = [openSet,neigcur(q)];
            openSet = orderize(openSet);
        elseif tent_gS >= gScore(neigcur(q))
            continue
        end
        cameFrom(neigcur(q)) = current;
        gScore(neigcur(q)) = tent_gS;
        %gScore
        fScore(neigcur(q)) = gScore(neigcur(q)) + satdistance(neigcur(q),s2,cfg);
        %fScore
    end
    
end

if ~exist('PATH', 'var')
    error('Oh shit!');
    %return;
end

L = length(PATH);
for q = 2:L
%     ind = [PATH(q-1),PATH(q)];
%     USAGE (min(ind),max(ind)) = 1;
    USAGE (PATH(q-1),PATH(q)) = 1;
end

end

