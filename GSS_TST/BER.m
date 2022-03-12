clear all;
R = 4000e3; %���������� ����� ��������� � ������
c = 299792458;
% f1=2.5e9;
f1=437e6;% �������
f2=24e9;
P1=3; %�������� ���������
P2=5;
df1=19.2e3; %������ ��������� �������
df2=56e3;
df3=1e6;


Antenna=4.5+4.5;%��������� �� �������� � ���������� ������
Prop=20*log10(c/f1)-20*log10(4*pi)-20*log10(R)+10*log10(2);% �������������� ����������� ���������� ()
%
Noise=-2-4+204-10*log10(df2);%���� ������������ ��������� �������� ����
L=Prop+Noise+Antenna; %����� �� ����������
Bqpsk=berawgn(L,'psk',4,'nondiff');%������� ����������� ������ � ����������� �� ���� ���������
Bfsk1=berawgn(L,'fsk',2,'noncoherent');
Bfsk2=berawgn(L,'fsk',2,'coherent');
