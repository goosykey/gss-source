function [ccoef, scoef] = readgm(fname)

% read gravity model data file

% input

%  fname = name of gravity data file

% output

%  ccoef, scoef = gravity model coefficients

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read the data file

gdata = dlmread(fname);

nrows = size(gdata, 1);

% initialize coefficients

idim = gdata(nrows, 1) + 1;

ccoef = zeros(idim, idim);

scoef = zeros(idim, idim);

% create gravity model coefficients

for n = 1:nrows
    
    i = gdata(n, 1);
    
    j = gdata(n, 2);
    
    ccoef(i + 1, j + 1) = gdata(n, 3);
    
    scoef(i + 1, j + 1) = gdata(n, 4);
    
end




