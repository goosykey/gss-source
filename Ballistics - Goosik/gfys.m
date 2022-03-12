singleframe(t1(1),[R1(:,1),R2(:,1),R3(:,1)]);
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');

boob  = 500000;

for k=1:10000 %length(t1)
    singleframe(t1(k+boob),[R1(:,k+boob),R2(:,k+boob),R3(:,k+boob)]);
    f=getframe;
    im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
end

imwrite(im,map,'imagefile.gif','DelayTime',0,'LoopCount',inf);