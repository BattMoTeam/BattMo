function coupterm = getCoupTerm(coupterms, coupname)
    
    counames = cellfun(@(c) c.name, coupterms, 'uniformoutput', false);
    ind = strcmp(coupname, coupnames);
    assert(any(ind), 'coupling term not found');
    assert(nnz(ind) == 1, 'several coupling terms with same name');
    coupterm = coupterms{ind(1)};
    
end