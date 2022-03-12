lats1 = 90:-1:-90;
lons1 = -180:179;
[latmap1, lonmap1] = ndgrid(lats1,lons1);
Omaskfull = oceanALL(latmap1,lonmap1,5);

weat1 = latmap1*0;
[weat2,weatb] = getweather(weat1,Omaskfull);
weatb = double(weatb);
weatb(~weatb)=nan;
surf(lons1,lats1,weatb)
view(0,90);
pause(0.5);

q = 30;
coco = zeros(1,q);



for i = 1:55
    [weat2,weatb] = getweather(weat2,Omaskfull);
    badperc = 1e2*sum(sum(weatb))/numel(weatb);
    goodperc = 1e2*sum(sum(~weatb))/numel(weatb);
    coco(i) = badperc;
    toptopstring = sprintf( 'BAD %3.2f - %3.2f GOOD',badperc,goodperc);
    weatb = double(weatb);
    cla reset
    weatb(~weatb)=nan;
    surf(lons1,lats1,weatb)
    text(0.5,1,toptopstring,'Units','normalized', 'VerticalAlignment','bottom',...
        'HorizontalAlignment','center', 'FontSize',10, 'Color','black');
    view(0,90);
    drawnow;
    pause(0.3);
end

%plot_google_map('maptype','satellite');

% figure
% plot(coco);

clear weat1 weat2 weatb goodperc badperc coco lats1 lons1 latmap1 lonmap1