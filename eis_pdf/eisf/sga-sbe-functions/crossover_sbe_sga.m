
% Two point crossover  

function [binpop,sga] = crossover_sbe_sga(binpop,sga)
    popsize = sga.popsize; 
    pcross_bin = sga.pcross_bin; 
    nbincross = sga.nbincross;
    bitlength = sga.bitlength; 
    for i = 1:2:popsize
        parent1 = binpop(i,:);
        parent2 = binpop(i+1,:);
        child1 = parent1;
        child2 = parent2;
        if rand < pcross_bin
            nbincross = nbincross + 1;
            point1 = randi([1 bitlength-2],1);
            point2 = randi([point1+1 bitlength-1],1);
            child1(point1:point2) = parent2(point1:point2);
            child2(point1:point2) = parent1(point1:point2);            
        end
        binpop(i,:) = child1; 
        binpop(i+1,:) = child2;        
    end    
    sga.nbincross = nbincross;
end
