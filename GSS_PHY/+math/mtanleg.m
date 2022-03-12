function [ mtanleg ] = mtanleg( lm,mm,fi )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mtP = zeros(lm+1,mm+1);

%mtP(:,1) = 0;


mtP(2,2) = sin(fi); %P11

for l = 2:lm
    for m = 1:l
        realleg = legendre(l-1,sin(fi)); %with kondon
        mtP(l+1,m+1) = mtP(l-1,m+1) + (2*l-1)*m*sin(fi)*(-1)^(m-1)*realleg(m);
    end
end

mtanleg = mtP;


end

