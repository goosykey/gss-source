%% INI

iniCDS;

throw_vel = 1e-3; % km/s

vel_error = 0.05; % error
alpha_error = 18*deg;

a0 = orb(1);
e0 = orb(2);
om0 = orb(3);
Om0 = orb(4);
in0 = orb(5);
u0 = orb(6);
m0 = orb(7);

[R0,V0] = math.randv(a0,e0,in0,Om0,om0,u0-om0);

VS = zeros(4,3);

for i = 1:4
    
    switch i
        case 1
            v_er = throw_vel * vel_error;
            alpha_er = 0;
        case 2
            v_er = -throw_vel * vel_error;
            alpha_er = 0;
        case 3
            v_er = throw_vel * vel_error;
            alpha_er = 0;
        otherwise
            v_er = -throw_vel * vel_error;
            alpha_er = 0;
    end
    
    
    v_er_vec = v_er * [cos(alpha_er)*cos(alpha_er) ...
        cos(alpha_er)*sin(alpha_er) sin(alpha_er)];
    v_er_vec = math.tsw2xyz(v_er_vec', Om0, in0, u0);
    
    VS(i,:) = V0 - throw_vel - v_er_vec;
end

[a1,e1,om1,Om1,in1,u1] = math.xyz2kepler(R0(1),R0(2),R0(3),VS(1,1),VS(1,2),VS(1,3));
[a2,e2,om2,Om2,in2,u2] = math.xyz2kepler(R0(1),R0(2),R0(3),VS(2,1),VS(2,2),VS(2,3));
[a3,e3,om3,Om3,in3,u3] = math.xyz2kepler(R0(1),R0(2),R0(3),VS(3,1),VS(3,2),VS(3,3));
[a4,e4,om4,Om4,in4,u4] = math.xyz2kepler(R0(1),R0(2),R0(3),VS(4,1),VS(4,2),VS(4,3));

Y1 = math.orbitize([a1,e1,om1,Om1,in1,u1,m0]);
Y2 = math.orbitize([a2,e2,om2,Om2,in2,u2,m0]);
Y3 = math.orbitize([a3,e3,om3,Om3,in3,u3,m0]);
Y4 = math.orbitize([a4,e4,om4,Om4,in4,u4,m0]);

%% test

%[ t,y,yfil ] = PEG( jd0, h, t_fin, dayperiod, Y, Y, 'A');

%supermegaplot(t, y, yfil, 'a');

%% BULL



for i = 1:4
    switch i
        case 1
            Y = Y1;
            Sbody = 0.09;
        case 2
            Y = Y2;
            Sbody = 0.09;
        case 3
            Y = Y3;
            Sbody = 0.01;
        case 4
            Y = Y4;
            Sbody = 0.01;
    end
    
    [ t,y,yfil ] = PEG( jd0, h, t_fin, dayperiod, Y, Y, 'A');
    
    y(diff(t)==0,:) = [];
    yfil(diff(t)==0,:) = [];
    t(diff(t)==0) = [];
    
    if i == 1
        ephcartfull = zeros(4,numel(t),6);
    end
    
    p = y(:,1);
    l1 = y(:,2);
    l2 = y(:,3);
    e = sqrt(l1.^2 + l2.^2);
    a = p./(1-e.^2);
    om = atan2(l2,l1);
    
    [Ri,Vi] = math.randv(a,e,y(:,5),y(:,4),om,y(:,6)-om);
    
    ephcartfull (i,:,:) = permute([Ri;Vi],[3,2,1]);
    
end

[R0,V0] = math.randv(a0,e0,in0,Om0,om0,u0-om0);

%% calc plot

ephc = ephcartfull(:,:,1:3);
ephcd = ephc*0;
ephcd (1,:,:) = ephc(2,:,:) - ephc(1,:,:);
ephcd (2,:,:) = ephc(3,:,:) - ephc(1,:,:);
ephcd (3,:,:) = ephc(4,:,:) - ephc(1,:,:);
ephcd (4,:,:) = ephc(3,:,:) - ephc(2,:,:);
ephcd (5,:,:) = ephc(4,:,:) - ephc(2,:,:);
ephcd (6,:,:) = ephc(4,:,:) - ephc(3,:,:);

dists = sqrt(sum(ephcd.^2,3));
distmax = max(dists,[],1);

plot (t/86400, distmax, '-b', 'linewidth', 2);
hold on; grid on;
plot (t/86400, distmax/2, '-g', 'linewidth', 2);

load('tomsk.mat');
plot (t/86400, distmaxtomsk, '--m', 'linewidth', 2);
plot (t/86400, distmaxtomsk/2.2, '--r', 'linewidth', 2);
legend('Container (10%)', 'Container (5%)', 'Manual deployment - rough', 'Manual deployment - precise');
xlabel('time, days');
ylabel('Separation distance, km');

%% VIS

movie_che(t, ephcartfull);










