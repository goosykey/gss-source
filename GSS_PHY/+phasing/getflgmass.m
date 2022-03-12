function [ flagmass ] = getflgmass( mainangdeg, polar )
%GET FLAG ARRAY
%   Returns an array of flags to use for thrust on/off
% 
%   INPUTS:
%       mainangdeg : Arg. of latitude in 1st quadrant to start or finish
%       polar      : boolean, 
%             if true  - thrust begins @mainangdeg degrees
%             if false - thrust ends @mainangdeg degrees
%  



flagmass = zeros(360,1);
us = (1:360)' ;

flagmass(us < mainangdeg | us > (360-mainangdeg)) = 1;   
flagmass(us > (180-mainangdeg) & us < (180+mainangdeg)) = 1;

flagmass = boolean(flagmass);

if polar
    flagmass = ~flagmass;
end

end

