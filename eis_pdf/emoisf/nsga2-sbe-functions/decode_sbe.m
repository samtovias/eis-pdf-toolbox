
% Decode for SBE 
% not decoding at all, just inserts the binary string in parents pop matrix

function pop = decode_sbe(pop, pop_strings, nobj, ncon, nreal)

    pop(:,nobj+ncon+nreal+1:end-3) = pop_strings;  

end
