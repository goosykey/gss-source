function plotem (hexasx, hexasy)

l = length(hexasx(:,1));

hexasx = [hexasx,hexasx(:,1)];
hexasy = [hexasy,hexasy(:,1)];

for i = 1:l
    plot(hexasx(i,:),hexasy(i,:)); hold on;
end

plot_google_map;

end