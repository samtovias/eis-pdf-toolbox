% Selection with binary tournament

function [bparents,id] = selection_gop(pop,bpop,params)
np = params.np;
bparents = false(size(bpop));
id = zeros(np,1);
% Generate random indices 
index = randperm(np);
% Add the first index to the final
index = [index index(1)];
% Example for 4 individuals:
% [3 4 1 2 3]: 3 vs 4, 4 vs 1, 1 vs 2, 2 vs 3 
for i=1:np
    find1 = pop(index(i),1);
    find2 = pop(index(i+1),1);
    % Binary tournament
    idx = tournament_gop(find1,find2);
    if idx == 1  
        bparents(i,:) = bpop(index(i),:);
        id(i) = index(i);
    else  
        bparents(i,:) = bpop(index(i+1),:);
        id(i) = index(i+1);
    end
end
end

% Binary tournament for evolutionary global optimization 
% function id = tournament_gop(find1,find2)
% if find1 < find2  
%     id = 1; 
% elseif find2 < find1  
%     id = 2; 
% elseif (rand <= 0.5)
%     id = 1; 
% else
%     id = 2; 
% end 
% end
