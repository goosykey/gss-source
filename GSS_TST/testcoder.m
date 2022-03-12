function [ali, asi] = testcoder(qcatz, shi_shoon)
%% this is to test cpp stuff

[~,N] = size(qcatz);
ja = numel(shi_shoon);

if N ~= ja
    disp('feck\n');
    ali = pi;
else
    ali = qcatz*shi_shoon';
end

asi = numel(ali);


end

