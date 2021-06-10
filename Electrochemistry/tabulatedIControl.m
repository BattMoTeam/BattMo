function [val, ctrltype] = tabulatedIControl(t, tpoints, Ipoints)

    ctrltype = 'I';
    val = interpTable(tpoints, Ipoints, t);
    
end
