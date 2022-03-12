%1 rot earth

jd0 = astro.jday(2017,3,20,10,28,0);

cdata = imread('EM_low.jpg');

%mins = 1440/4;

for i = 0:5:1440
    jd = jd0+i/1440;
    cla;
    vis_earthdraw(jd,'axes','both','cdata',cdata);
    view(0,45);
    pause(0.01);
end