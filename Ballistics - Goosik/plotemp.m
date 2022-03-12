figure('Name','Report: Semi-major axis','NumberTitle','off');
plot(t1/86400, eph1(:,1), 'r', 'LineWidth', 2); grid; hold;
plot(t2/86400, eph2(:,1), 'g', 'LineWidth', 2);
plot(t3/86400, eph3(:,1), 'b', 'LineWidth', 2);
legend('a_1(t)','a_2(t)','a_3(t)', 0);
xlabel('time, days');
ylabel('Semi-major axis, km, deg');

figure('Name','Report: RAAN','NumberTitle','off');
plot(t1/86400, eph1(:,4)*180/pi, 'r', 'LineWidth', 2); grid; hold;
plot(t2/86400, eph2(:,4)*180/pi, 'g', 'LineWidth', 2);
plot(t3/86400, eph3(:,4)*180/pi, 'b', 'LineWidth', 2);
legend('\Omega_1(t)','\Omega_2(t)','\Omega_3(t)', 0);
xlabel('time, days');
ylabel('RAAN, deg');

figure('Name','Report: INCLINATION','NumberTitle','off');
plot(t1/86400, eph1(:,5)*180/pi, 'r', 'LineWidth', 2); grid; hold;
plot(t2/86400, eph2(:,5)*180/pi, 'g', 'LineWidth', 2);
plot(t3/86400, eph3(:,5)*180/pi, 'b', 'LineWidth', 2); hold;
legend('i_1(t)','i_2(t)','i_3(t)', 0);
xlabel('time, days');
ylabel('Inclination, deg');

figure('Name','Report: MASS','NumberTitle','off');
plot(t1/86400, eph1(:,7), 'r', 'LineWidth', 2); grid; hold;
plot(t2/86400, eph2(:,7), 'g', 'LineWidth', 2);
plot(t3/86400, eph3(:,7), 'b', 'LineWidth', 2);
legend('m_1(t)','m_2(t)','m_3(t)', 0);
xlabel('time, days');
ylabel('Mass, kg');