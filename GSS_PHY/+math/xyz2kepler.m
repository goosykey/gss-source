function [ a,e,om,Om,in,u ] = xyz2kepler( x,y,z,vx,vy,vz )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = 398600.4418; % �������������� �������� ����� (= G*M), ��3/�2 <- ��������� ��������, ��������� ���
rxyz = [x,y,z]'; % ������ ������
vxyz = [vx,vy,vz]'; % ������ �������� ��

Ixyz = cross(rxyz,vxyz); % �����-�� ��� ������ ��������
I = norm(Ixyz); % ��� ������

r = norm([x,y,z]); % ������ ������-�������
v = norm([vx,vy,vz]); % ������ ������� ��������

Esp = v^2/2 - mu/r; % ������������� �������

a = -mu/(2*Esp); % ������� ������� ������
e = sqrt(1-I^2/(a*mu)); % �������������� ������
p = a*(1-e^2); % ��������� �������� ������

in = acos(Ixyz(3)/I); % ����������

Om = atan2(Ixyz(1),-Ixyz(2)); % ������� ����������� ����

u = atan2(z/sin(in), x*cos(Om)+y*sin(Om)); % �������� ������

if abs(p-r) > 1E-04 % ���� ������ �� ��������
    nu = atan2(sqrt(p/mu)*dot(vxyz,rxyz),p-r); % �������� �������� �������� ��
    om = u - nu; % �������� ���������� ������
else
    nu = 0; % ����� ��� ��� �������� �� ����� ����������� ������,
    om = 0; % �� ��� ���� �� �������� ���������� �������
end

end

