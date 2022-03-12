function [ penta_c, pentas, hexa_c, hexas, F_unused, allF ] = hexxing( F )
%HEXXING - build a penta-hex mesh on a triangular mesh
%   Detailed explanation goes here

Flast = max(max(F));

penta_c = [];
pentas = [];
hexa_c = [];
hexas = [];

used_triangles = [];
centres_overall = [];

%% FIND PENTAGONS

fprintf('Finding pentagons\n');

counts = sum(hist(F,max(max(F))),2); % how many times a vertice is in F
pentacentres = [counts,(1:numel(counts))'];
pentacentres = pentacentres(pentacentres(:,1)==5,2);

for i = pentacentres'
    Bhere = any(F==i,2);
    used_triangles = [used_triangles;F(Bhere,:)];
    %number = sum(Bhere);
    
    %fprintf('PENTA!! %3.2f > ',i);
    penta_c = [penta_c ; i];
    Fhere = F(Bhere,:);
    
    vec = Fhere(1,:);
    vec = vec(vec~=i); vec = vec(end);
    Fhere(1,:)=[];
    for j = 1:4
        lastv = vec(end);
        newno = any(Fhere==lastv,2);
        newel = Fhere(newno,:);
        newel = newel(newel~=lastv & newel ~= i);
        vec = [vec,newel];
        Fhere(newno,:)=[];
    end
    
    pentas = [pentas;vec];
    F(Bhere,:) = nan;
    
end


centres_overall = pentacentres;

%% HEXAGONS

fprintf('Finding hexagons\n');
% prompt = 'Go? >';
% command = input (prompt);

counter = 0;

while numel(F(~isnan(F)))/3 > 0
    try
        counter = counter + 1;
        fprintf('Hexagon generation... step %g  || %g triangles left\n',counter,numel(F(~isnan(F)))/3);
        tic
        centersofused = ismember (used_triangles, centres_overall);
        triangles_modif = used_triangles;
        triangles_modif(centersofused) = 0;
        triangles_modif = sort(triangles_modif,2);
        used_pairs = triangles_modif(:,2:3);
        
        counts = sum(hist(F,max(max(F))),2); % how many times a vertice is in F
        validpoints = [counts,(1:numel(counts))'];
        validpoints = validpoints(ismember(validpoints(:,1),[1,2,3,4]),2);
        
        ism_pairs = ismember (used_pairs,validpoints);
        ism_pairs = sum(ism_pairs,2);
        
        validpairs = used_pairs(ism_pairs==2,:);
        fprintf('Valid points : %g ; Valid pairs : %g\n',numel(validpoints),numel(validpairs)/2);
        
        ism_1 = ismember(F,validpairs(:,1));
        ism_2 = ismember(F,validpairs(:,2));
        ism_0 = ism_1 | ism_2;
        numbers = sum(ism_0,2) == 2;
        
        FF = F(numbers,:);
        ism_1 = ismember(FF,validpairs(:,1));
        ism_2 = ismember(FF,validpairs(:,2));
        ism_0 = ~(ism_1 | ism_2);
        hexacentres_raw = FF(ism_0);
        hexacentres = hexacentres_raw(~ismember(hexacentres_raw,centres_overall));
        hexacentres = orderize(hexacentres);
        
        if isempty(hexacentres)
            break
        end
        
        fprintf('Hexacentres found : %g ; Unused of them : %g\n',numel(hexacentres_raw),numel(hexacentres));
        
        for i = hexacentres'
            Bhere = any(F==i,2);
            used_triangles = [used_triangles;F(Bhere,:)];
            %number = sum(Bhere);
            
            %fprintf('PENTA!! %3.2f > ',i);
            hexa_c = [hexa_c ; i];
            Fhere = F(Bhere,:);
            
            if isempty(Fhere)
                continue
            end
            
            vec = Fhere(1,:);
            vec = vec(vec~=i); vec = vec(end);
            Fhere(1,:)=[];
            for j = 1:5
                lastv = vec(end);
                newno = any(Fhere==lastv,2);
                newel = Fhere(newno,:);
                newel = newel(newel~=lastv & newel ~= i);
                vec = [vec,newel];
                Fhere(newno,:)=[];
            end
            
            if numel(vec) < 6
                NN = 6 - numel(vec);
                for k = 1:NN
                    vec = [vec,nan];
                end
            end
            
            hexas = [hexas;vec];
            F(Bhere,:) = nan;
            
        end
        
        centres_overall = [pentacentres;hexacentres];
        
        toto = toc;
        fprintf('Triangles left: %g ; Time taken this step: %3.1f sec.\n',numel(F(~isnan(F)))/3, toto);
        F = F(~isnan(F(:,1)),:);
        fprintf('\n');
        
    catch
        vec
        break
    end
    
    %     prompt = 'Go on? > ';
    %     command = input (prompt);
    
end

F_unused = F;

pentasq = zeros(length(pentas(:,1)),1);
pentasq(:) = nan;



end

