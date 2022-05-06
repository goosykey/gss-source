function [ N ] = Nlm( l, m )
%NNM Summary of this function goes here
%   Normalization of gravity constants and Legendre functions
%   Pn = N*P (Legendre); C = N*Cn (Constants)

del = (m==0);

N = sqrt( factorial(l-m) * (2*l+1) * (2 - del)/factorial(l+m) );


end