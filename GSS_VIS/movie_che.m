function fig = movie_che(t, ephcart, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   ephcart EXPECTED CUBE (Nsat, Nt, 6)

%% CONST
deg = pi/180;
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

timeres = 0.02; % THIS IS REAL TIME RESOLUTION

%% PRE-DEF

timescale = 240; % X REAL TIME
jd0 = 2.4599455e+06;
satstyles = 'or';
satcolors = [0.75 0.75 0.75];

orbitpersistence = 43200;

%% VARARGIN CHECK

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'jd0'
            jd0 = varargin{2};
        case 'satstyles'
            satstyles = varargin{2};
        case 'satcolors'
            satcolors = varargin{2};
        case 'timescale'
            timescale = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%%DEFINITIONS & CALCULATIONS

[Nsat, ~, ~] = size(ephcart);

dt = t(2) - t(1);

dtnew = timeres * timescale; % model time step while 0.04 real seconds pass
tnew = t(1):dtnew:t(end);
Ntnew = numel(tnew);

ephcartnew = zeros(Nsat,numel(tnew),3);

for i = 1:3
    ephcartnew(:,:,i) = interp1(t,ephcart(:,:,i)',tnew)'; % NEW EPH DATA
end

globerotangle = erot*dtnew/deg; % IN DEGREES

%% IMPL

[ fig, globe, misc ] = vis_earthdraw(jd0, 'quality','low');

X = ephcartnew(:,1,1); Y = ephcartnew(:,1,2); Z = ephcartnew(:,1,3);
sats = plot3(X,Y,Z,satstyles, 'markerfacecolor', satcolors);

XL1 = X([1 2]); YL1 = Y([1 2]); ZL1 = Z([1 2]); L1 = plot3(XL1,YL1,ZL1,'-c', 'linewidth', 3);
XL2 = X([1 3]); YL2 = Y([1 3]); ZL2 = Z([1 3]); L2 = plot3(XL2,YL2,ZL2,'-c', 'linewidth', 3);
XL3 = X([1 4]); YL3 = Y([1 4]); ZL3 = Z([1 4]); L3 = plot3(XL3,YL3,ZL3,'-c', 'linewidth', 3);
XL4 = X([2 3]); YL4 = Y([2 3]); ZL4 = Z([2 3]); L4 = plot3(XL4,YL4,ZL4,'-c', 'linewidth', 3);
XL5 = X([2 4]); YL5 = Y([2 4]); ZL5 = Z([2 4]); L5 = plot3(XL5,YL5,ZL5,'-c', 'linewidth', 3);
XL6 = X([3 4]); YL6 = Y([3 4]); ZL6 = Z([3 4]); L6 = plot3(XL6,YL6,ZL6,'-c', 'linewidth', 3);

toplstring = sprintf('t = + %3.2f day \n%3.2f\n%3.2f\n%3.2f\n%3.2f\n%3.2f\n%3.2f\n',...
    t(1)/86400,...
    sqrt( (X(2)-X(1)).^2 + (Y(2)-Y(1)).^2 + (Z(2)-Z(1)).^2) , ...
    sqrt( (X(3)-X(1)).^2 + (Y(3)-Y(1)).^2 + (Z(3)-Z(1)).^2) , ...
    sqrt( (X(4)-X(1)).^2 + (Y(4)-Y(1)).^2 + (Z(4)-Z(1)).^2) , ...
    sqrt( (X(3)-X(2)).^2 + (Y(3)-Y(2)).^2 + (Z(3)-Z(2)).^2) , ...
    sqrt( (X(4)-X(2)).^2 + (Y(4)-Y(2)).^2 + (Z(4)-Z(2)).^2) , ...
    sqrt( (X(4)-X(3)).^2 + (Y(4)-Y(3)).^2 + (Z(4)-Z(3)).^2) );

perssteps = round(orbitpersistence/dtnew);

X = ephcartnew(:,1:2,1); Y = ephcartnew(:,1:2,2); Z = ephcartnew(:,1:2,3);
orbits = plot3(X',Y',Z','-r');


% tuxt = text(0,1,toplstring,'Units','normalized', 'VerticalAlignment','top',...
%     'HorizontalAlignment','left', 'FontSize',10, 'Color','red');

an = annotation(fig, 'textbox', [0 0 0 0], 'String', toplstring, 'FitBoxToText','on',  'verticalalignment', 'bottom', 'color',[1 1 1]);

for i = round(Ntnew/2):Ntnew
    rotate([globe; misc],[0 0 1], globerotangle);
    X = ephcartnew(:,i,1); Y = ephcartnew(:,i,2); Z = ephcartnew(:,i,3);
    sats.XData = X;
    sats.YData = Y;
    sats.ZData = Z;
    
    %view(atan2(Y(1),X(1))/deg,asin(Z(1)/sqrt(X(1)^2+Y(1)^2+Z(1)^2))/deg)
    %view([X(1) Y(1) Z(1)]);
    
    XL1 = X([1 2]); YL1 = Y([1 2]); ZL1 = Z([1 2]);
    XL2 = X([1 3]); YL2 = Y([1 3]); ZL2 = Z([1 3]);
    XL3 = X([1 4]); YL3 = Y([1 4]); ZL3 = Z([1 4]);
    XL4 = X([2 3]); YL4 = Y([2 3]); ZL4 = Z([2 3]);
    XL5 = X([2 4]); YL5 = Y([2 4]); ZL5 = Z([2 4]);
    XL6 = X([3 4]); YL6 = Y([3 4]); ZL6 = Z([3 4]);
    
    
    L1.XData = XL1; L1.YData = YL1; L1.ZData = ZL1;
    L2.XData = XL2; L2.YData = YL2; L2.ZData = ZL2;
    L3.XData = XL3; L3.YData = YL3; L3.ZData = ZL3;
    L4.XData = XL4; L4.YData = YL4; L4.ZData = ZL4;
    L5.XData = XL5; L5.YData = YL5; L5.ZData = ZL5;
    L6.XData = XL6; L6.YData = YL6; L6.ZData = ZL6;
    
    toplstring = sprintf('t = +%3.2f day\n1 > 2 : %3.2f\n1 > 3 : %3.2f\n1 > 4 : %3.2f\n2 > 3 : %3.2f\n2 > 4 : %3.2f\n3 > 4 : %3.2f\n',...
        tnew(i)/86400,...
        sqrt( (X(2)-X(1)).^2 + (Y(2)-Y(1)).^2 + (Z(2)-Z(1)).^2) , ...
        sqrt( (X(3)-X(1)).^2 + (Y(3)-Y(1)).^2 + (Z(3)-Z(1)).^2) , ...
        sqrt( (X(4)-X(1)).^2 + (Y(4)-Y(1)).^2 + (Z(4)-Z(1)).^2) , ...
        sqrt( (X(3)-X(2)).^2 + (Y(3)-Y(2)).^2 + (Z(3)-Z(2)).^2) , ...
        sqrt( (X(4)-X(2)).^2 + (Y(4)-Y(2)).^2 + (Z(4)-Z(2)).^2) , ...
        sqrt( (X(4)-X(3)).^2 + (Y(4)-Y(3)).^2 + (Z(4)-Z(3)).^2) );
    
    %tuxt.String = toplstring;
    an.String = toplstring;
    
    if mod(i,1000) == 0 || i == 2 || i == round(Ntnew/2)  % Plot orbits
        if i + perssteps/2 <= Ntnew
            if i - perssteps/2 >= 1
                persrange = round(i - perssteps/2) : round (i + perssteps/2);
            else
                persrange = 1 : round (i + perssteps/2);
            end
        else
            if i - perssteps/2 >= 1
                persrange = round(i - perssteps/2) : Ntnew;
            else
                persrange = 1 : Ntnew;
            end
        end
    end
    X = ephcartnew(:,persrange,1); Y = ephcartnew(:,persrange,2); Z = ephcartnew(:,persrange,3);
    delete(orbits);
    orbits = plot3(X',Y',Z','-r');
    pause (timeres);
    
end


end

