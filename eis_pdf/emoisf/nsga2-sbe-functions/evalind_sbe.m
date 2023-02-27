
% Evaluate individual with SBE (nbits instead nbin)

function pop = evalind_sbe(pop, popsize, nobj, ncon, nreal, nbits,...
                                     X, Y, params) 

xbin = pop(:, nobj+ncon+nreal+1:nobj+ncon+nreal+nbits); 

if (nreal == 0) && (nbits == 0) 
    disp('There are no binary and real variables defined!');
    return;
end
 
% Instance selection with simple binary encoding 
f1 = zeros(popsize,1);
f2 = zeros(popsize,1); 
parfor i=1:popsize
    [f1(i),f2(i)] = evalind_sbe_capsule(X,Y,xbin(i,:),params);
end
pop(:,1) = f1; 
pop(:,2) = f2; 

for i=1:popsize
    if (ncon == 0)
        pop(:,nobj+ncon+nreal+nbits+1) = 0.0;
        break
    else
        pop(i,nobj+ncon+nreal+nbits+1) = 0.0;
        for j=1:ncon
            if (pop(i,nobj+j) < 0.0)
                pop(i,nobj+ncon+nreal+nbits+1) = pop(i,nobj+ncon+nreal+nbits+1) + pop(i,nobj+j);
            end
        end
    end
end
