aa1(:,:) = aa (end:-1:1,:);

s = size(aa1);

CC = {[],[],[],[],[],[],[]};
figure; hold on;

for k = 1:7
    for i = 1:s(1)
        for j = 1:s(2)
            if aa1(i,j) == k
                CC{k} = [CC{k}; [i j]];
            end
        end
    end
end

plot(CC{1}(:,2),CC{1}(:,1),'+r');
plot(CC{2}(:,2),CC{2}(:,1),'+b');
plot(CC{3}(:,2),CC{3}(:,1),'xm');
plot(CC{4}(:,2),CC{4}(:,1),'sg');
plot(CC{5}(:,2),CC{5}(:,1),'ok','markerfacecolor','r');
plot(CC{6}(:,2),CC{6}(:,1),'*c');
plot(CC{7}(:,2),CC{7}(:,1),'.y');

legend('Monopropellant (H_4N_2)', 'Monopropellant (H_2O_2)',...
    'Bipropellant (NTO/MMH)', 'Arcjet (NH3)',...
    'Ion (Xe)', 'VAT (Cr)', 'FEEP');