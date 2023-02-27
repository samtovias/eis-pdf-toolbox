
% Selection with binary tournament

function newbinpop = selection_sga(binpop,fitpop)
    popsize = size(binpop,1);
    ind = single((1:popsize)');
    % Binary tournament
    vs = [Shuffle(ind) Shuffle(ind)];     
    winners = fitpop(vs(:,1)) <= fitpop(vs(:,2)); 
    index = [vs(winners,1);vs(~winners,2)];
    newbinpop = binpop(index,:);
end
