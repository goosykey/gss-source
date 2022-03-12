function [B, eps] = access2(r, rg, eps0)
%ACCESS1 Summary of this function goes here
%   Detailed explanation goes here

r1 = r - rg;

r1 = r1 ./ repmat(sqrt(sum(r1.^2,2)),[1,3]);

rg = rg ./ repmat(sqrt(sum(rg.^2,2)),[1,3]);

eps = asind(dot(rg,r1,2));

B = eps > eps0;



end

