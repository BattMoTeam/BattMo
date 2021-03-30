function coupterm = getCoupTerm(coupterms, coupname, coupnames)
    
    ind = strcmp(coupname, coupnames);
    assert(any(ind), 'coupling term not found');
    assert(nnz(ind) == 1, 'several coupling terms with same name');
    coupterm = coupterms{ind};
    
end