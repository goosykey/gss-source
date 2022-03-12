rr = 3.0;
freq = [1:1000]*1e9;
L = rainpl(10000,freq,rr);
loglog(freq/1e9,L)
grid
xlabel('Frequency (GHz)')
ylabel('Attenuation (dB)')

