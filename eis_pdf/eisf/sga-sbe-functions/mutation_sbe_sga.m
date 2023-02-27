
% Mutation  

function [binpop,sga] = mutation_sbe_sga(binpop,sga)
    pmut_bin = sga.pmut_bin;
    nbinmut = sga.nbinmut; 
    if pmut_bin == -1; pmut_bin = 1/sga.bitlength; end
    idx = rand(size(binpop)) < pmut_bin;
    nbinmut = nbinmut + sum(sum(idx));
    sga.nbinmut = nbinmut;
    binpop(idx) = ~binpop(idx);
end 