
% Decode binary population 

function realpop = decode_sga(binpop,sga)
popsize = sga.popsize; 
min_binvar = sga.min_binvar; 
max_binvar = sga.max_binvar;
nbits = sga.nbits;
nbin = sga.nbin;
realpop = zeros(popsize,numel(nbits));
for i=1:popsize
    stop = 0;
    for j=1:nbin
        start = stop+1; stop = start + nbits(j) - 1;
        total = sum(binpop(i,start:stop) .* 2.^(nbits(j)-1 : -1: 0));
        realpop(i,j) = min_binvar(j) + (max_binvar(j)-min_binvar(j)) / (2^nbits(j)-1) * total;
    end
end