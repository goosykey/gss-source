function plotCircle3D(center,normal,radius, formatstring)

if nargin < 4
    formatstring = '-r';
end

theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),formatstring);

end