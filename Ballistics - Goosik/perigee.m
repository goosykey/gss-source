function [ om ] = perigee( l1, l2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if l1 == 0
    om = pi/2 * sign(l2);
elseif l1 > 0
    om = atan(l2./l1);
else
    om = atan(l2./l1)+pi;
end

om(om(:)<0) = om(om(:)<0) + 2*pi;
    

end

