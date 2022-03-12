function [ output_args ] = plotusage( SATCOOR, CONN, cfg, USAGE, DEPOTS, chhubs )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin<4
    plotbrokenlinks = false;
end

pl = cfg(1); spp = cfg(2);
satcou = pl * spp;
SATCOOR1 = SATCOOR;
SATCOOR1(spp:spp:satcou,2) = SATCOOR1(spp:spp:satcou,2)-360;

CONN_1 = real(CONN).*heaviside(real(CONN)) + imag(CONN).*heaviside(imag(CONN))*1i;
CONN_0 = real(CONN).*heaviside(-real(CONN)) + imag(CONN).*heaviside(-imag(CONN))*1i;

global temp temp1 temp2;
temp = CONN;
temp1 = CONN_1;
temp2 = CONN_0;


usage = USAGE;

USAGE = USAGE + USAGE.';

figure('Name','Plot','NumberTitle','on'); hold on;

levels = round(max(max(USAGE)));
if levels <= 1
    huestep = 1;
else
    huestep = 240/(levels-1)/360;
end
paletteHSV = zeros(levels,3);
for q = 1:levels
    paletteHSV(q,:) = abs([2/3-(q-1)*huestep, 1, 1]);
end

%paletteHSV
paletteRGB = hsv2rgb(paletteHSV);


plot(SATCOOR(:,1),SATCOOR(:,2), 'o');hold on;

XLim = get(gca,'XLim'); XX = XLim(2)-XLim(1);
YLim = get(gca,'YLim'); YY = YLim(2)-YLim(1);
pos = get(gcf,'Position'); PP = pos(4)/pos(3);

for q = 1:satcou
    for w = 1:satcou
        if w <= q
            continue
        end
        
        if CONN_1(q,w) == 0
            continue;
        end
        
        if USAGE(q,w) ==0
            if real(CONN_1(q,w)) > 0
                plot([SATCOOR(q,1),SATCOOR(w,1)],[SATCOOR(q,2),SATCOOR(w,2)],'k','LineWidth',1);
            else
                plot([SATCOOR1(q,1),SATCOOR1(w,1)],[SATCOOR1(q,2),SATCOOR1(w,2)],'k','LineWidth',1);
            end
            continue
        end
        
        stringmode = (usage(w,q)~=0)-(usage(q,w)~=0);
        switch stringmode
            case 0
                %stri = ['\leftarrow',num2str(usage(w,q)),' || ',num2str(usage(q,w)),' \rightarrow'];
                stri = sprintf('\\leftarrow %3.2f || %3.2f \\rightarrow',usage(w,q),usage(q,w));
            case 1
                %stri = ['\leftarrow',num2str(usage(w,q))];
                stri = sprintf('\\leftarrow %3.2f',usage(w,q));
            case -1
                %stri = [num2str(usage(q,w)),' \rightarrow'];
                stri = sprintf('%3.2f \\rightarrow',usage(q,w));
        end
        
        if real(CONN_1(q,w)) ~= 0
            x1 = SATCOOR(q,1);
            x2 = SATCOOR(w,1);
            y1 = SATCOOR(q,2);
            y2 = SATCOOR(w,2);
            rotangle = atan((y2-y1)/(x2-x1)*XX/YY*PP)/pi*180;
            if round(USAGE(q,w)) == 0
                plot([x1,x2],[y1,y2],'Color',paletteRGB(round(USAGE(q,w)+1),:),'LineWidth',2);
            else
                plot([x1,x2],[y1,y2],'Color',paletteRGB(round(USAGE(q,w)  ),:),'LineWidth',2);
            end
            text((x1+x2)/2, (y1+y2)/2,...
                stri,...
                'FontSize', 8, 'HorizontalAlignment', 'center', 'Rotation', rotangle);
            continue
        end
        
        if imag(CONN_1(q,w)) ~= 0
            x1 = SATCOOR1(q,1);
            x2 = SATCOOR1(w,1);
            y1 = SATCOOR1(q,2);
            y2 = SATCOOR1(w,2);
            rotangle = atan((y2-y1)/(x2-x1)*XX/YY*PP)/pi*180;
            if round(USAGE(q,w)) == 0
                plot([x1,x2],[y1,y2],'Color',paletteRGB(round(USAGE(q,w)+1),:),'LineWidth',2);
            else
                plot([x1,x2],[y1,y2],'Color',paletteRGB(round(USAGE(q,w)  ),:),'LineWidth',2);
            end
            text((x1+x2)/2, (y1+y2)/2,...
                stri,...
                'FontSize', 8, 'HorizontalAlignment', 'center', 'Rotation', rotangle);
        end
    end
end

% gplot(real(CONN_1),SATCOOR, '-g');
% gplot(imag(CONN_1),SATCOOR1,'-g');

gplot(real(CONN_0),SATCOOR, ':k');
gplot(imag(CONN_0),SATCOOR1,':k');

for q = 1:length(DEPOTS)
    plot(SATCOOR(DEPOTS(q),1),SATCOOR(DEPOTS(q),2), 'mp', 'MarkerSize', 16, 'MarkerFaceColor', 'm');
end

for q = 1:satcou
    if any (DEPOTS==q)
        text(SATCOOR(q,1),SATCOOR(q,2), num2str(q),'HorizontalAlignment', 'center', 'FontSize',8);
    else
        text(SATCOOR(q,1),SATCOOR(q,2), [num2str(q),'>',num2str(chhubs(q))],...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize',8);
    end
end

%axis square;




end

