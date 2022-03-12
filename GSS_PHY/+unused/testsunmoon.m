function [ MSM ] = testsunmoon( N )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


muM = 7.36E+22 * 6.67384E-20; % Moon grav.param. km3/s2
muS = 1.98892E+30 * 6.67384E-20; % Sun grav.param. km3/s2
au1 = 149597871; % 1 Astronomic Unit in kms
R0 = 6378.14; % mean Earth radius

a = 3000;
b = 6000;

reci = zeros(N,3);
fsunmoonxyz = zeros(N,3);
jd = zeros(N,1);



for i = 1:N
    clc
    fprintf('%3.2f%% COMPLETE;\n', (i-1)/N*100);
    sig = sign(rand-0.5);
    reci(i,:) = (a + (b-a)*rand(1,3))*sig;
    jd(i) = 2451545.0 + 12000 * rand;
    if i == 1
        reci(i,:) = [5780.57872015653, 3024.88697217792, 5473.88502799918];
        jd(i) = 2460753.03041633;        
    end
    
    [rsun,rtascsun,declsun] = sun ( jd(i) ); %Vallado
    
    r12 = moon(jd(i))'*R0; r12abs = sqrt(sum(r12 .* r12));
    r13 = rsun'*au1; r13abs = sqrt(sum(r13 .* r13));
    r2 = reci(i,:)' - r12;    r2abs = sqrt(sum(r2 .* r2));
    r3 = reci(i,:)' - r13;    r3abs = sqrt(sum(r3 .* r3));
    
    %XYZ cartesian 3body perturbations
    fsunmoonxyz(i,:) = -muM/r2abs^3*r2 - muM/r12abs^3*r12 ...
        - muS/r3abs^3*r3 - muS/r13abs^3*r13;
end

MSM = [reci,jd, fsunmoonxyz];

end

