function [ a_fin, in_fin ] = RGTSSO( norbits, ndays, ini_inc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    ini_inc = 95;
end

TOL = 0.001; %km
R0 = 6378.137; %km

iter = 2;

[a, in] = RGTfixinc(norbits, ndays, ini_inc);
[a_new, in_new] = RGTfixinc(norbits, ndays, in);


while abs(a_new-a) > TOL
    iter = iter+1;
    a = a_new;
    in = in_new;
    [a_new, in_new] = RGTfixinc(norbits, ndays, in);
end

a_fin = a_new;
in_fin = in_new;

fprintf('======= CONVERGED =======   \n\n');

fprintf('number of iterations               %12.0f   \n', iter);

fprintf('mean semimajor axis fin            %12.6f  kilometers \n', a_fin);

fprintf('mean altitude fin                  %12.6f  kilometers \n', a_fin-R0);

fprintf('mean inclination fin               %12.6f  degrees \n\n', in_fin);


end

