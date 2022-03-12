function [ VS, DS, SS, Pout, VSd, DSd, SSd, VSa, DSa, SSa] = getener_old( a, e, om, Om, in, u, netdata, SESDATA, SESCOORD, jd )
%function [ NS, DS, SS, Pout, NSd, DSd, SSd ] = getener( a, e, om, Om, in, u, SESDATA, SESCOORD, utc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isempty(SESDATA)
    VS = 0; DS = 0; SS = 0;
    VSa = 0; DSa = 0; SSa = 0;
    VSd = 0; DSd = 0; SSd = 0;
    Pout = 0;
    return
end

dOm = netdata{1};
du = netdata{2};
ddu = netdata{3};

deg = pi/180;
R0 = 6378.14;
sector = 50*deg;
h = a - R0;

c = 3e5; % speed of light km/s
f_down = 20e9; %Hz
f_up = 17e9; % Hz

losses_km = 0.07;
losses_km_up = 0.04;

SNR_m = 27;
SNR_f = 27;
Gain_AT_m = 12;
Gain_AT_f = 23;

k = 1.381*1e-23;
T = 500;
N0 = 10*log10(k*T);

Losses_tech = -7.9;

C_R=1./2;
K_CP=1+1./4;

bitrate_m = netdata{10};
bitrate_f = netdata{11};

% vtau = netdata{4};
% dtau = netdata{5};
% stau = netdata{6};
% 
% utcsec = ceil(mod(jd-0.5,1) * 86400);
%utcsec = 36005;

% SESUTCSEC = ceil(mod(SESBEG-0.5,1) * 86400);
% 
% INTOSESSION = utcsec-SESUTCSEC;
% 
% SESUTCSEC(INTOSESSION(:) < 0) = NaN; % sessions not begun yet
% SESUTCSEC(boolean((INTOSESSION(:) > vtau) .* (SESDATA(:,1)==1))) = NaN; % v.sessions already over
% SESUTCSEC(boolean((INTOSESSION(:) > dtau) .* (SESDATA(:,1)==2))) = NaN; % d.sessions already over
% SESUTCSEC(boolean((INTOSESSION(:) > stau) .* (SESDATA(:,1)==3))) = NaN; % s.sessions already over
% 
% 
% SESDATA(isnan(SESUTCSEC(:)),:) = [];
% SESCOORD(isnan(SESUTCSEC(:)),:) = [];
% SESUTCSEC(isnan(SESUTCSEC(:))) = [];

A = [a a a a a a a];
E = [e e e e e e e];
PERIG = [om om om om om om om];
OMEGA = [Om Om Om+dOm Om+dOm Om Om-dOm Om-dOm];
IN = [in in in in in in in];
U = [u u+ddu u+du u+du-ddu u-ddu u-du u-du+ddu];

[Rxyz,~] = randv(A,E,IN,OMEGA,PERIG,U-PERIG);

Sxyz = [R0.*cosd(SESCOORD(:,1)').*cosd(SESCOORD(:,2)');R0.*cosd(SESCOORD(:,1)').*sind(SESCOORD(:,2)');R0.*sind(SESCOORD(:,1)')];

P_all = zeros(1,7);
Sass = [0,0,0];

for i = 7:-1:1
    Ri = Rxyz(:,i);
    Ri = R3(gstime(jd))*Ri;
    %ri = norm(Ri);
    %fprintf('lat %g lon %g ||',asin(Ri(3)/ri)/deg,atan2(Ri(2),Ri(1))/deg);
    Ritab = zeros(1,numel(SESCOORD)/2);
    for ii = 1:3
        Ritab(ii,:) =  Ri(ii);
    end
    DISTANCES = Sxyz;
    RiRi = Ritab;
    
    DISTANCES = (DISTANCES - RiRi).^2;
    
    DISTANCES = sum(DISTANCES);
    
    DISTANCES = sqrt(DISTANCES);
    
    SiRi = Sxyz - Ritab;
    
    SiRi_normalizer = [sqrt(sum(SiRi.^2));sqrt(sum(SiRi.^2));sqrt(sum(SiRi.^2))];
    SiRi1 = SiRi./SiRi_normalizer;
        
    THETAS = acos(SiRi1'*(-Ri)/(1.0*norm(Ri)));
    THETAS(THETAS(:) > sector) = nan;
    DISTANCES(isnan(THETAS(:))) = nan;
    DISTANCES(DISTANCES(:) > 3000) = nan;
    THETAS(isnan(DISTANCES(:))) = nan;
    %DISTANCES(isnan(SESUTCSEC(:))) = nan;
    %THETAS(isnan(SESUTCSEC(:))) = nan;
    %SESUTCSEC(isnan(THETAS(:))) = NaN;
    
    SESDATA(~isnan(DISTANCES(:)),2) = SESDATA(~isnan(DISTANCES(:)),2) + 2^(i-1); % ADD SATELLITE FLAG TO SESDATA
    DISTANCES = DISTANCES'; % TRANSPOSE ROWv -> COLv
    
    
    fs_losses = ((4*pi*DISTANCES*f_down/c).^2);
    %expansion = fs_losses/min(fs_losses(~isnan(fs_losses(:))));
    fs_losses = -10*log10(fs_losses);
    
    alpha_introp = asin(((h+R0)./(R0+20)).*sin(THETAS))-THETAS;
    distance_introp = DISTANCES - (R0+20).*sin(alpha_introp)./sin(THETAS);
    distance_introp(isnan(distance_introp(:))) = 20;
    Losses_introp = -losses_km * distance_introp;
    
    Gain0 = 50.9;
    Pol_losses = 10*log10(cos(THETAS));
    Gain = Gain0+0.3.*Pol_losses+10.*log10(cos(THETAS));
    
    PSD_m = -Gain_AT_m - Gain - fs_losses - Losses_introp - Losses_tech + N0 + SNR_m;
    PSD_d = -Gain_AT_m - Gain - fs_losses - Losses_introp - Losses_tech + N0 + SNR_m;
    PSD_f = -Gain_AT_f - Gain - fs_losses - Losses_introp - Losses_tech + N0 + SNR_f;
    
    Ptx_m = PSD_m + 10*log10(bitrate_m*K_CP/C_R);
    Ptx_m(SESDATA(:,1)~=1) = nan;
    
    Ptx_d = PSD_d + 10*log10(bitrate_m*K_CP/C_R);
    Ptx_d(SESDATA(:,1)~=2) = nan;
    
    Ptx_f = PSD_f + 10*log10(bitrate_f/4*(1+1/16)/(7/8));
    Ptx_f(SESDATA(:,1)~=3) = nan;
    
    
    Ptx_m_W = 10.^(0.1.*Ptx_m);
    Ptx_d_W = 10.^(0.1.*Ptx_d);
    Ptx_f_W = 10.^(0.1.*Ptx_f);
    
    PV = sum(Ptx_m_W(~isnan(Ptx_m_W)));
    PD = sum(Ptx_d_W(~isnan(Ptx_d_W)));
    PS = sum(Ptx_f_W(~isnan(Ptx_f_W)));
    P_all(i) = PV+PD+PS;
    
    if i > 1
        continue
    end
    
    VSd = sum(~isnan(Ptx_m_W));
    DSd = sum(~isnan(Ptx_d_W));
    SSd = sum(~isnan(Ptx_f_W));
    
    POWAHFULL = [(1:numel(Ptx_m_W))',Ptx_m_W,Ptx_d_W,Ptx_f_W,double(SESDATA(:,2))];
    POWAHFULL(mod(POWAHFULL(:,5),2) == 0,:) = []; % sessions reachable by FIRST and any other satellite(s) or 1st only
    POWAH = POWAHFULL(POWAHFULL(:,5) > 1,:); % sessions reachable by FIRST and any other satellite(s)
    %global temp temp1 temp2 temp3;
    %temp1 = []; temp2 = []; temp3 = [];
    %temp = {POWAHFULL,POWAH};
    
    while P_all(1) > 700 && ~isempty(POWAH); % 'handover' sequence
        POWAHMAX = max(POWAH(:,2:4),[],2);
        [maxp,maxi] = max(POWAHMAX);
        %temp1 = [temp1;maxp];
        indexh = POWAH(maxi,1);
        %temp2 = [temp2;maxi,indexh];
        %temp3 = [temp3;P_all];
        typeh = SESDATA(indexh,1);
        flagh = de2bi(POWAH(maxi,5),7);
        assfound = false;
        for q = 2:7
            if flagh(q) && P_all(q) < 700 % assist found
                assfound = true;
                %fprintf('%g',q);
                P_all(q) = P_all(q) + maxp;
                P_all(1) = P_all(1) - maxp;
                Sass(typeh) = Sass(typeh) + 1;
                switch typeh % session handed over to a neighbour
                    case 1
                        Ptx_m_W(indexh) = nan;
                    case 2
                        Ptx_d_W(indexh) = nan;
                    case 3
                        Ptx_f_W(indexh) = nan;
                end
%                 POWAHFULL(indexh,:)=[];
                POWAH(maxi,:)=[];
                break % stop looking for assist
            end
        end  
        if ~assfound % handover loop failed, assist not found for this sub
            %fprintf('N');
%             POWAHFULL(POWAHFULL(:,1)==indexh,5) = 1; % session scheduled for possible drop
%             POWAH = POWAHFULL(POWAHFULL(:,5) > 1,:);
            POWAH(maxi,:) = [];
        end
    end
    
    while P_all(1) > 700 && sum(~isnan(Ptx_d_W)) > 0 % drop dataN subs
         %fprintf('D');
         
         [maxp,maxid] = max(Ptx_d_W);
         Ptx_d_W(maxid) = nan;
         
         P_all(1) = P_all(1) - maxp;
     end
    
    while P_all(1) > 700 % drop maximum power
         %fprintf('D');
         [maxpv,maxiv] = max(Ptx_m_W);
         [maxpd,maxid] = max(Ptx_d_W);
         [maxps,maxis] = max(Ptx_f_W);
         
         [maxp,maxi] = max([maxpv,maxpd,maxps]);
         
         switch maxi
             case 1
                 Ptx_m_W(maxiv) = nan;
             case 2
                 Ptx_d_W(maxid) = nan;
             case 3
                 Ptx_f_W(maxis) = nan;
         end
         P_all(1) = P_all(1) - maxp;
     end
    
    
    VS = sum(~isnan(Ptx_m_W)); VSa = Sass(1) + VS;
    DS = sum(~isnan(Ptx_d_W)); DSa = Sass(2) + DS;
    SS = sum(~isnan(Ptx_f_W)); SSa = Sass(3) + SS;
    Pout = P_all(1);
     
end



end

