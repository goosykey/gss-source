function [ Nreal, Pout, Ntot, Ndes, Nass] = getener( a, e, om, Om, in, u, netdata, SESDATA, SESCOORD, jd )
%function [ NS, DS, SS, Pout, NSd, DSd, SSd ] = getener( a, e, om, Om, in, u, SESDATA, SESCOORD, utc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Ntot =  [VSt, DSt, SSt];
% Ndes = [VSd, DSd, SSd];
% Nass = [VSa, DSa, SSa];
% Nreal = [VS, DS, SS];

if isempty(SESDATA)
    Ntot =  [0, 0, 0];
    Ndes = [0, 0, 0];
    Nass = [0, 0, 0];
    Nreal = [0, 0, 0];
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

C_R=1./2; % coderate
K_CP=1+1./4; % cyclic preffix

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
Sdro = [0,0,0];

sesq = length(SESDATA(:,1));
P_full = zeros(sesq,11);
P_full(:,10) = 1:sesq;

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
    
    Ptx_W = [Ptx_m_W(SESDATA(:,1) ==  1) ; Ptx_d_W(SESDATA(:,1) ==  2) ; Ptx_f_W(SESDATA(:,1) ==  3) ];
    P_full(:,i) = Ptx_W;
    
    if i > 1
        continue
    end
    
    P_full(:,8:9) = double(SESDATA);
    
    global temp temp1 temp2 temp3;
    temp = P_full;
    temp3 = SESDATA;
    
    %     numel(~isnan(Ptx_m_W))
    %     numel(~isnan(Ptx_d_W))
    %     numel(~isnan(Ptx_f_W))
    
    VSt = sum(~isnan(Ptx_m_W));
    DSt = sum(~isnan(Ptx_d_W));
    SSt = sum(~isnan(Ptx_f_W));
    Ntot =  [VSt, DSt, SSt];
    
    P_full(P_full(:,9)==0,:) = []; % delete sessions not belonging any sat coverage zone
    
    if isempty(P_full)
        Ntot =  [0, 0, 0];
        Ndes = [0, 0, 0];
        Nass = [0, 0, 0];
        Nreal = [0, 0, 0];
        Pout = 0;
        return
    end
    
    sesq_new = length(P_full(:,1));
    P_full(:,11) = 1:sesq_new; % add numbers in cropped P_full table
    
    [P_chosen, SAT_chosen] = min(P_full(:,1:7),[],2);
    
    
    
    temp1 = P_chosen;
    temp2 = SAT_chosen;
    flag_chosen = 2.^(SAT_chosen-1);
    flag_chosen_bi = de2bi(flag_chosen,7);
    
    P_notfull = P_full;
    P_notfull(:,1:7) = P_notfull(:,1:7) .* flag_chosen_bi;
    P_notfull(P_notfull == 0) = nan;
    P_notfull(:,9) = flag_chosen;
    
    VSd = sum(~isnan(P_notfull(P_notfull(:,8)==1,1)));
    DSd = sum(~isnan(P_notfull(P_notfull(:,8)==2,1)));
    SSd = sum(~isnan(P_notfull(P_notfull(:,8)==3,1)));
    Ndes = [VSd, DSd, SSd]; % DESIRED SUBS after first distribution
    
    P_all = sum(P_notfull(:,1:7),1,'omitnan'); % preliminary P_out distribution estimation
    
    P_assistable = P_full(P_notfull(:,9)==1 & P_full(:,9) > 1,:); % sessions picked by sat#1 able to pass to any other sat
    
    %     ASSISTABLE = sum(P_notfull(:,9)==1 & P_full(:,9) > 1 ,1);
    ASSISTABLE = length(P_assistable(:,1));
    
    
    
    while P_all(1) > 700 && ASSISTABLE > 0
        [phere,seshere] = max(P_assistable(:,1),[],1);
        seshere_full = P_assistable(seshere,11); % number of the session in P_full table
        seshere_true = P_assistable(seshere,10); % true number of the session
        flaghere = de2bi(P_assistable(seshere,9),7);
        typehere = P_assistable(seshere,8); % type of the session
        assfound = false;
        for q = 2:7
            if flaghere(q) && P_all(q) < 700 % assist found
                assfound = true;
                %fprintf('%g',q);
                P_all(q) = P_all(q) + P_full(seshere_full,q);
                P_all(1) = P_all(1) - phere;
                
                P_notfull(seshere_full,1) = nan;
                P_notfull(seshere_full,q) = P_full(seshere_full,q);
                P_notfull(seshere_full,9) = 2^(q-1);
                
                P_assistable(seshere,:) = [];
                ASSISTABLE = ASSISTABLE-1;
                
                Sass(typehere) = Sass(typehere) + 1;
                
                break % stop looking for assist
            end
        end
        if ~assfound % assist loop failed, this is not an assistable session
            %fprintf('N');
            P_assistable(seshere,:) = []; % only delete this session from the assistables table
            ASSISTABLE = ASSISTABLE - 1;
        end
    end
    
    P_notfull(mod(P_notfull(:,9),2) == 0,:) = []; % get rid of sessions served by other sats
    P_notfull(:,12) = 1:length(P_notfull(:,1)); % numbers in new P_notfull table
    
    P_droppable = P_notfull(P_notfull(:,8)==2 & P_notfull(:,9)==1 , :); % extract only dataN sessions from sat#1
    DROPPABLE = length(P_droppable(:,1));
    
    while P_all(1) > 700 && DROPPABLE > 0 % drop dataN subs
        %fprintf('D');
        
        [phere,seshere] = max(P_droppable(:,1));
        seshere_notfull = P_droppable(seshere,12); % number of the session in P_notfull table
        
        P_all(1) = P_all(1) - phere;
        
        P_notfull(seshere_notfull,1) = nan;
        
        P_droppable(seshere,:) = [];
        DROPPABLE = DROPPABLE-1;
        
        Sdro(2) = Sdro(2) + 1;
        
    end
    
    P_droppable = P_notfull(P_notfull(:,8)~=2 & P_notfull(:,9)==1 , :); % extract all other sessions from sat#1
    
    while P_all(1) > 700 % drop maximum power
        %fprintf('D');
        [phere,seshere] = max(P_droppable(:,1));
        seshere_notfull = P_droppable(seshere,12); % number of the session in P_notfull table
        typehere = P_droppable(seshere,8);
        
        P_all(1) = P_all(1) - phere;
        
        P_notfull(seshere_notfull,1) = nan;
        P_droppable(seshere,:) = [];
        Sdro(typehere) = Sdro(typehere) + 1;
    end
    
    
    VS = VSd - Sass(1) - Sdro(1); VSa = Sass(1) + VS;
    DS = DSd - Sass(2) - Sdro(2); DSa = Sass(2) + DS;
    SS = SSd - Sass(3) - Sdro(3); SSa = Sass(3) + SS;
    Pout = P_all(1);
    Nass = [VSa, DSa, SSa];
    Nreal = [VS, DS, SS];
    
end



end

