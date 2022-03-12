function [ rmagn ] = magnitude( Rxyz )
%Returns magnitudes of vectors
%   Returns 1xN array of magnitudes of a MxN input column-vector array
%   (typically M = 3)

L = length(Rxyz(1,:));

rmagn = zeros (1,L);

for i = 1:L
    rmagn(i) = norm(Rxyz(:,i));
end


end

