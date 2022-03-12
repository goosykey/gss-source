function [delv,epsv, deltas, deltaa ] = navdelta( a,e,om,Om,in,u, delmax, epsmax, maxind )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[R0,V0] = randv(a,e,in,Om,om,u-om);

delv = -delmax:delmax/(maxind-1):delmax;
epsv = -epsmax:epsmax/(maxind-1):epsmax;


p = a*(1-e^2);

L = length(delv);
deltaa = zeros(L);
deltas = zeros(L^2,8);

for i = 1 : L
    for j = 1 : L
        
        R = R0+delv(i);
        V = V0+epsv(j);
        
        x=R(1); y=R(2); z=R(3);
        vx=V(1); vy=V(2); vz=V(3);
        
        [ a1,e1,om1,Om1,in1,u1 ] = xyz2kepler2( x,y,z,vx,vy,vz );
        
        p1 = a1*(1-e1^2);
        
        da = a1-a;
        de = e1-e;
        dom = om1-om;
        dOm = Om1-Om;
        din = in1-in;
        du = u1-u;
        
        dp = p1-p;
        
        if dom > pi && dom < 2*pi
            dom = dom-2*pi;
        end
        
        deltas(L*(i-1)+j,:) = ([delv(i)*1e3,epsv(j)*1e3, da*1e3,de,dom,dOm,din,du ]);
        deltaa(i,j) = da*1E03;
    end
end



end

