function plotcone( CART510, alpha, colour )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(CART510(:,1));


for i = 1:N
    add3d.elem.covcone( CART510(i,1:3), alpha, colour )
end


end

