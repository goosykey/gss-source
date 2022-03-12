function [ output_args ] = plotrelv( t, rxyz, vxyz, fi, la )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[ rrelgnd, vrelgnd, vrelgndrad, visible ] = relgnd( t, rxyz, vxyz, fi, la );

rrelmag = magnitude(rrelgnd);
vrelmag = magnitude(vrelgnd);

for i = 1:length(rrelmag)
    if ~(visible(i)==1)
        vrelmag(i) = nan;
        vrelgndrad(i) = nan;
    end
end

figure('Name','Relative Velocity RADIAL','NumberTitle','off');
plot (t/86400, abs(vrelgndrad), 'r', 'LineWidth', 2);
grid;
xlabel('time, days');
ylabel('relative velocity, m/s');

figure('Name','Visibility','NumberTitle','off');
plot (t/86400, visible, 'g', 'LineWidth', 4);
grid;
xlabel('time, days');
ylabel('VISIBLE');



end

