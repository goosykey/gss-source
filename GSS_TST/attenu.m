freq = [1:200]*1e9;

R = 1000.0;

T = 15;
P = 101300.0;
W = 7.5;
L = gaspl(R,freq,T,P,W);

L0 = gaspl(R,freq,T,P,0.0);

semilogy(freq/1e9,L)
hold on
semilogy(freq/1e9,L0)
grid
xlabel('Frequency (GHz)')
ylabel('Specific Attenuation (dB/km)')
hold off

legend('Air + H_2O','Dry air');