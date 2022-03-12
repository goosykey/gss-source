function [ mprop ] = tsiolk( m0, deltav, Ispsec )
%TSIOLK Summary of this function goes here
%   Detailed explanation goes here

mprop = m0 * (exp(deltav/Ispsec/9.81)-1);

end

