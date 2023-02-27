
% Binary tournament for evolutionary global optimization 
function id = tournament_gop(find1,find2)
if find1 < find2  
    id = 1; 
elseif find2 < find1  
    id = 2; 
elseif (rand <= 0.5)
    id = 1; 
else
    id = 2; 
end 
end

% function selected = selection_gop(bpop,fpop)
% 
%     np = size(bpop,1);
%     ind = (1:np)';
%     vs = [Shuffle(ind) Shuffle(ind)];  
%     while sum(vs(:,1)==vs(:,2)) > 0
%         vs = [Shuffle(ind) Shuffle(ind)];
%     end
%     winners = fpop(vs(:,1)) <= fpop(vs(:,2)); 
%     index = [vs(winners,1);vs(~winners,2)];
%     selected = bpop(index,:);
% 
% end

% j1 = 0; j2 = 0;
% while j1==j2 % Ensure different parents
%     j1 = TournamentSelection(Fpop,params);
%     j2 = TournamentSelection(Fpop,params);
% end
