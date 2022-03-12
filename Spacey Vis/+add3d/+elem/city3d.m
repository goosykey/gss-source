function city3d( jd, lat, lon, names, varargin )
%CITY3D Summary of this function goes here
%   VECTOR INPUT SUPPORTED

color = 'm';
markersize = 4;
fontsize = 7;


while ~isempty(varargin)
    switch lower(varargin{1})
        case 'markersize'
            markersize = varargin{2};
        case 'fontsize'
            fontsize = varargin{2};
        case 'color'
            color = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

R0 = 6378.137;

thetag = astro.gstime(jd);

recf = [R0.*cosd(lat).*cosd(lon) R0.*cosd(lat).*sind(lon) R0.*sind(lat)]';

reci = math.R3(-thetag)*recf;

reci = reci * 1e3 * (R0+50)/R0;

plot3(reci(1,:),reci(2,:),reci(3,:),'om', 'markerfacecolor', color, 'markersize',markersize, 'markeredgecolor', 'none');

text(reci(1,:),reci(2,:),reci(3,:),names,'VerticalAlignment','bottom',...
        'HorizontalAlignment','center', 'color',color, 'FontSize',fontsize);


end

