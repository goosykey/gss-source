function [ val_out ] = looptrapz( qq, t1, t2 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if t1 == t2
    val_out = trapz(qq);
elseif t1 < t2
    val_out = trapz(qq(t1:t2));
else
    val_out = trapz(qq)-trapz(qq(t2:t1));
end


end

