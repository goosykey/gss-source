function [ SATCOOR, CONNout ] = generatenet( cfg, faultprob, aol )
%Model ISL traffic

deg = pi/180;
aol = aol/deg;
aol(aol>180) = aol(aol>180) - 360;
aol = abs(aol);
aol(aol>90) = aol(aol>90) - 180;
aol = abs(aol);

pl = cfg(1); spp = cfg(2);
satcou = pl * spp;

%% LIMITATIONS

aolimit = 75; %AOL LIMIT

BRIDGENO = [0, round(spp/4), round(2*spp/4), round(3*spp/4)];

%%

if nargin < 3 || isempty(aol)
    aol = zeros(satcou,1);
end






SATCOOR = zeros(satcou,2);

for q = 1:satcou
    planeno = idivide(q-1,int32(spp))+1;
    satinplane = mod(q,spp);
    if satinplane == 0
        satinplane = spp;
    end
    SATCOOR(q,:) = [planeno, (planeno-1)*(360/spp*2/3) + (satinplane-1)*360/spp];
end

SATCOOR1 = SATCOOR;
SATCOOR1(spp:spp:satcou,2) = SATCOOR1(spp:spp:satcou,2)-360;

CONN = zeros(satcou);

for q = 1:satcou
    for w = 1:satcou
        if (w-q)==1 && mod(q,spp)~=0 %FORWARD
            CONN(q,w)=1;
        end
        
        if (w-q)==spp && aol(q) < aolimit && any(BRIDGENO==mod(q,spp))
            CONN(q,w)=1;
        end
        
%         if (w-q)==spp-1 && aol(w) < aolimit && any(q,BRIDGENO)
%             CONN(q,w)=1;
%         end
%         
%         if (w-q)==2*spp-1 && mod(w,spp)==0 && aol(q) < aolimit && any(q,BRIDGENO)
%             CONN(q,w)=1i;
%         end
        
        if (w-q)==spp-1 && mod(w,spp)==0 %FORWARD IM
            CONN(q,w)=1i;
        end
    end
end

for q = 1:satcou
    for w = 1:satcou
        if rand <= faultprob
            CONN(q,w)=CONN(q,w)*-1;
        end
    end
end
CONNout = CONN;

end

