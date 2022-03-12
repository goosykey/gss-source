function plotcov( CART510, alpha, formatstring )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(CART510(:,1));

for i = 1:N
    add3d.elem.covcircle( CART510(i,1:3), alpha, formatstring )
end


end

