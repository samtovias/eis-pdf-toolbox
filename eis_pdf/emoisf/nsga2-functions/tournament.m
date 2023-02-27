function parent = tournament(ind1, ind2, nobj, ncon, nreal, nbin)
flag = check_dominance(ind1, ind2, nobj, ncon, nreal, nbin);
if (flag == 1)
    parent = ind1;
elseif (flag == -1)
    parent = ind2;
elseif (ind1(1,end) > ind2(1,end))
    parent = ind1;
elseif (ind2(1,end) > ind1(1,end))
    parent = ind2;
elseif (rand <= 0.5)
    parent = ind1;
else
    parent = ind2;
end
