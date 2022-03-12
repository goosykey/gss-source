clear
clc

lat = -90:5:90;
tt = zeros(1,numel(lat));

[~, ~, ~, ~,~, ~, ephcart, gscart] = ...
            access1('time', 0:10:20e04); % initialize ephemeris
        

for i = 1:numel(lat)
    fprintf('|');
    [~, ~, tt(i), ~] = access1('poi',[lat(i),0], 'time', 0:10:20e04, 'error', 0, 'ephcart', ephcart);
end

fprintf('\n');

figure;
plot(lat,tt/60);
hold on;
grid;

tt = tt'/60;