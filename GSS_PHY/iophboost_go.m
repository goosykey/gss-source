ga = 40:140;
massess = zeros(numel(ga),1);
timess = zeros(numel(ga),1);
tete = cell(numel(ga),1);

m_disp = 550;

figure;

for i = 1:numel(ga)
    N = ga(i);
    [ masses, times ] = iophboost(N,m_disp,0);
    massess(i) = masses(1);
    timess(i) = times(end);
    tete{i} = sprintf('N = %g',N);
end

massess = massess - 980 - 100 - m_disp - 4200;

plot(timess/3600,massess,'*-', 'linewidth',2);
hold on;
grid on;
xlabel('time, hours');
ylabel('propellant mass, kg');
text(timess/3600,massess,tete,'VerticalAlignment','bottom');

plot(timess/3600,ones(1,numel(ga))*(600-m_disp),'-r', 'linewidth',2);

clear times masses timess massess ga N m_disp tete
