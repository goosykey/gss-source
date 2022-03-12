function [  ] = worldpop( map, BIGMASK )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[ Vmask, Dmask, Smask, ~, ~, Pmask, SPmask, ~, Zmask] = decompose( BIGMASK );

Nsubmap = map.*Pmask;
Ssubmap = map.*SPmask;

[M,N] = size(map);

Vmask = imresize(double(Vmask)/100,[M N],'bilinear');
Dmask = imresize(double(Dmask)/100,[M N],'bilinear');
Smask = imresize(double(Smask)/100,[M N],'bilinear');

Vmap = map.*Vmask.*Pmask;
Dmap = map.*Dmask.*Pmask;
Smap = map.*Smask.*SPmask;

Nsubscr = sum(sum(Nsubmap));
Ssubscr = sum(sum(Ssubmap));

Vmodsubs = sum(sum(Vmap));
Dmodsubs = sum(sum(Dmap));
Smodsubs = sum(sum(Smap));

VmodsubsZ = [sum(Vmap(Zmask(:)==0)) sum(Vmap(Zmask(:)==1)) sum(Vmap(Zmask(:)==2)) sum(Vmap(Zmask(:) >2))];

people = sum(sum(map));
people0 = round(sum(map(Zmask(:)==0))); Nsubscr0 = round(sum(Nsubmap(Zmask(:)==0))); Ssubscr0 = round(sum(Ssubmap(Zmask(:)==0)));
people1 = round(sum(map(Zmask(:)==1))); Nsubscr1 = round(sum(Nsubmap(Zmask(:)==1))); Ssubscr1 = round(sum(Ssubmap(Zmask(:)==1)));
people2 = round(sum(map(Zmask(:)==2))); Nsubscr2 = round(sum(Nsubmap(Zmask(:)==2))); Ssubscr2 = round(sum(Ssubmap(Zmask(:)==2)));
people3 = round(sum(map(Zmask(:) >2))); Nsubscr3 = round(sum(Nsubmap(Zmask(:) >2))); Ssubscr3 = round(sum(Ssubmap(Zmask(:) >2)));

figure('Name','POPULATION ZONAL DISTRIBUTION','NumberTitle','on'); hold on;
vec = [people0, people1, people2, people3];
labels = {'Distant','Rural','Urban','Ocean'};
pie (vec,labels);
hold off;axis off

figure('Name','SUBSCRIBER ZONAL DISTRIBUTION','NumberTitle','on'); hold on;

Z = [Nsubscr0+Ssubscr0, Nsubscr1+Ssubscr1, Nsubscr2+Ssubscr2, Nsubscr3+Ssubscr3];
labels = {'Distant','Rural','Urban','Ocean'};
ax1 = subplot(1,1,1);
pie(ax1,Z,labels)
title(ax1,'All subs');
hold off;axis off

figure('Name','VOICE USAGE ZONAL DISTRIBUTION','NumberTitle','on'); hold on;
vec = VmodsubsZ;
labels = {'Distant','Rural','Urban','Ocean'};
pie (vec,labels);
hold off;axis off

fprintf('In area: TOTAL (DISTANT + RURAL + URBAN + OCEAN) \n\tPopulation : %g thsnd (%g + %g + %g + %g)\n\tNorm subs  : %g (%g + %g + %g + %g)\n\tSpeed subs : %g (%g + %g + %g + %g)\n',...
    people/1000, people0/1000, people1/1000, people2/1000, people3/1000, Nsubscr, Nsubscr0, Nsubscr1, Nsubscr2, Nsubscr3, Ssubscr, Ssubscr0, Ssubscr1, Ssubscr2, Ssubscr3);
fprintf('Average intensity \n\tVoice : %3.2f\n\tDataN : %3.2f\n\tDataS : %3.2f\n', Vmodsubs/Nsubscr, Dmodsubs/Nsubscr, Smodsubs/Ssubscr);


end

