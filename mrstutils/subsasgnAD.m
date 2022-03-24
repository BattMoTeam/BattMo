function u = subsasgnAD(u,ind,v)
    if ~isa(u, 'ADI') && isa(v,'ADI')% u is a vector
            u = double2ADI(u, v);
    end
    u(ind) = v;             
end