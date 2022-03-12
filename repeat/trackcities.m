%% Moscow orbit - find cities
% recf expected: 3xn

function citylist = trackcities(jd0, t, recf, varargin)

deg = pi/180;
R0 = 6378.14;

% default
census = 1e6;
tolerance = 20; % km
litonly = 1;

nvarargs = length(varargin);

if rem(nvarargs,2) ~= 0
    error('Invalid number of arguments' )
end

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'census'
            census = varargin{2};
        case 'tolerance'
            tolerance = varargin{2};
        case 'litonly'
            litonly = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

cdb = readcities('cities1000.txt', census);        % read city database

cdb(:,6:7) = {[]};
cdb(:,8:12) = {[]};

tnew = 0:t(end); tnew = tnew';

[t1, index] = unique(t);

recf = recf';
recf = interp1(t1, recf(index,:), tnew,'spline');

lateph = asind(recf(:,3)./sqrt(sum(recf.^2,2)));
loneph = atan2d(recf(:,2),recf(:,1));

ncities = length(cdb(:,1));

%% CYCLE - loop = number of cities

fprintf('Cities - loop... \n');

for i = 1:ncities
    DISTS = deg*R0*distance(cdb{i,4}, cdb{i,5}, lateph, loneph);
    ifclose = DISTS < tolerance ;
    if ~any(ifclose)
        cdb(i,8) = {nan};
        continue
    end
    deltaifclose = ifclose(2:end) - ifclose(1:end-1);
    deltaifclose = [ifclose(1);deltaifclose]; % if begins with 1
    if ifclose(end) % if ends with 1
        deltaifclose(end) = -1;
    end
    cdb(i,6) = {tnew((deltaifclose==1))}; % points of commencement
    cdb(i,7) = {tnew((deltaifclose==-1))}; % points of termination
    cdb(i,8) = { numel( cdb{i,7} ) }; % number of passes
    cdb(i,9) = {cdb{i,7} - cdb{i,6}};
    cdb(i,10) = { sum(cdb{i,9}) };
    cdb(i,11) = { min(DISTS) };
    [~,day,~] = astro.gdate(jd0+tnew((deltaifclose==1))/86400);
    UTCangle = (day - floor(day))*360;
    LTangle = UTCangle + cdb{i,5};
    LTangle = mod(LTangle,360);
    cdb(i,12) = {LTangle/360*24};
end

jo = cdb(:,8);
jo = cell2mat(jo);
cdb(isnan(jo),:) = [];

if litonly
    ncities = length(cdb(:,1));
    for i = 1:ncities
        if numel(cdb{i,12}) == 1
            if cdb{i,12} < 7 || cdb{i,12} > 17
                cdb(i,12) = {nan};
            end
        else
            cdb{i,6}(cdb{i,12} < 7 | cdb{i,12} > 17) = [];
            cdb{i,7}(cdb{i,12} < 7 | cdb{i,12} > 17) = [];
            cdb{i,9}(cdb{i,12} < 7 | cdb{i,12} > 17) = [];
            cdb{i,12}(cdb{i,12} < 7 | cdb{i,12} > 17) = [];
        end
    end
    
    jo = cdb(:,12);
    jo = cell2mat(jo);
    cdb(isnan(jo),:) = [];
end


citylist = cdb;


end