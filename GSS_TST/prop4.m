clear; clc;

Isp = 100:0.1:500;
Mdry = 100;

mu = 398600.44;
R = 6878.14;
deg = pi/180;
g0 = 10;

v0 = sqrt(mu/R);

dv_man = v0*1e3*4*deg

Mp_sk = Mdry * (exp(100./Isp./g0)-1);

Mp_man = (Mdry+Mp_sk) .* (exp(dv_man./Isp./g0)-1);

M_total = 3*Mdry + 3*Mp_sk + 2*Mp_man;

%plot(Isp, M_total, 'r', 'linewidth',2); hold on; grid;

[~, mini] = min(abs(M_total-350));

I = Isp(mini)
M_man = Mp_man(mini)
M_sk = Mp_sk(mini)

