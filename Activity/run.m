%% GETTRAFFIC EXAMPLES

HR = 11.7; % HOURS OF DAY
TC = 10; % SECONDS
flg = 23; % PIE CHART FLAG
EX = 1; % TEST NUMBER SEE BELOW

%% SWITCH EXAMPLE

switch EX
    case 1 % NEW YORK, USA
        gettraffic( 40.709757, -74.009411, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 2 % AUCKLAND, NEW ZEALAND
        gettraffic( -36.897561, 174.762829, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 3 % MOSCOW, RUSSIA
        gettraffic( 55.750026, 37.624041, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 4 % TOKYO, JAPAN
        gettraffic( 35.662568, 139.770193, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 5 % COLORADO, USA, [RECTANGULAR EXAMPLE]
        gettraffic( [41.000402 36.993142], [-109.049929 -102.042039], 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );        
    case 6 % AUSTIN-DALLAS-HOUSTON, TX, USA
        gettraffic( 31.142518, -96.597478, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 7 % BELO HORIZONTE, BRAZIL
        gettraffic( -19.923536, -43.947772, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 8 % CORAL SEA : VANUATU & SOLOMON ISLES
        gettraffic( -12.775276, 164.512592, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 9 % SCANDINAVIA
        gettraffic( 62.553800, 19.836737, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 10 % SOUTH AFRICA
        gettraffic( -27.765062, 24.502797, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );   
    case 11 % NORTHERN SOUTH AMERICA
        gettraffic( 5.238286, -69.508834, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg ); 
    case 12 % WESTERN AFRICA
        gettraffic( 10.746479, -6.246217, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );    
    case 13 % NORTHERN AFRICA
        gettraffic( 31.632168, 1.918312, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );   
    case 14 % CENTRAL AFRICA
        gettraffic( -2.560044, 29.776781, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg ); 
    case 15 % SOUTHERN INDIA
        gettraffic( 15.719318, 76.341709, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
    case 16 % CENTRAL KAZAKHSTAN
        gettraffic( 48.927762, 70.141174, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );  
    case 17 % WESTERN EUROPE
        gettraffic( 48.355766, 7.961011, 750, HR, HR+TC/3600, map+map_oc, mapdata, BIGASSMASK, 0, flg );
end

clear EX HR TC flg;