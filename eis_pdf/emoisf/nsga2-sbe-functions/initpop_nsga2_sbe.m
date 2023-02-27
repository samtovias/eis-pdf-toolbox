
% Init binary population in NSGA2 with SBE 

function parent_strings = initpop_nsga2_sbe(popsize,bitlength,typeinit)

    parent_strings = zeros(popsize,bitlength);
    
    if strcmp('sparse',typeinit)
        
        min = 0.05; 
        max = 0.95;
        for i = 1:popsize
            u = (max-min)*rand + min;
            for j = 1:bitlength
                if u < rand
                    parent_strings(i,j) = 0;
                else 
                    parent_strings(i,j) = 1; 
                end
            end
        end
        
    else
        
        parent_strings = round(rand(popsize,bitlength));
    
    end

end
