function [K] = loveNumbers(n,m)
%LOVENUMBERS Summary of this function goes here
%   Detailed explanation goes here

switch n
    case 2
        switch m
            case 0
                K = [0.29525 -0.00087 0.30190 -0.00000 -0.00089];
            case 1
                K = [0.29470 -0.00079 0.29830 -0.00144 -0.00080];
            case 2
                K = [0.29801 -0.00057 0.30102 -0.00130 -0.00057];
            otherwise
                error('Wrong order! must be 0..2');
        end
    case 3
        switch m
            case 0
                K = [0.093 nan 0.093 0 nan];
            case 1
                K = [0.093 nan 0.093 0 nan];
            case 2
                K = [0.093 nan 0.093 0 nan];
            case 3
                K = [0.094 nan 0.094 0 nan];
            otherwise
                error('Wrong order! must be 0..3');
        end
    otherwise 
        error('Wrong degree! Must be 2 or 3');
end

end