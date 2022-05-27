classdef CompGraphFilter
    
    properties
        FieldToAddList
        FieldToRemoveList
    end
    methods
        
        function cgf = CompGraphFilter()
            
        end
        
        function cgf = addFieldToAdd(cgf, str)
            if isempty(cgf.FieldToAddList)
                cgf.FieldToAddList = {str};
            else
                cgf.FieldToAddList = {cgf.FieldToAddList{:}, str};
            end
        end
        
        function cgf = addFieldToRemove(cgf, str)
            if isempty(cgf.FieldToRemoveList)
                cgf.FieldToRemoveList = {str};
            else
                cgf.FieldToRemoveList = {cgf.FieldToRemoveList{:}, str};
            end
        end        
        
        
        function [varnames, propfuncs] = applyFilter(cgf, varnames, propfuncs)

            fdas = cgf.FieldToAddList;
            fdrs = cgf.FieldToRemoveList;
            
            if ~isempty(fdas)

                newvarnames = {};
                for ivarname = 1 : numel(varnames)
                    name = varnames{ivarname}.name;
                    for ifda = 1 : numel(fdas)
                        if regexp(name, fdas{ifda})
                            newvarnames{end + 1} = varnames{ivarname};
                            break
                        end
                    end
                end
                varnames = newvarnames;
                
                newpropfuncs = {};
                for iprop = 1 : numel(propfuncs)
                    name = propfuncs{iprop}.varname.name;
                    for ifda = 1 : numel(fdas)
                        if regexp(name, fdas{ifda})
                            newpropfuncs{end + 1} = propfuncs{iprop};
                            break
                        end
                    end
                end
                propfuncs = newpropfuncs;
                
            end
            
            
            if ~isempty(fdrs)

                newvarnames = {};
                for ivarname = 1 : numel(varnames)
                    name = varnames{ivarname}.name;
                    tokeep = true;
                    for ifdr = 1 : numel(fdrs)
                        if regexp(name, fdrs{ifdr})
                            tokeep = false;
                            break
                        end
                    end
                    if tokeep
                        newvarnames{end + 1} = varnames{ivarname};
                    end
                end
                varnames = newvarnames;

                newpropfuncs = {};
                for iprop = 1 : numel(propfuncs)
                    name = propfuncs{iprop}.varname.name;
                    tokeep = true;
                    for ifdr = 1 : numel(fdrs)
                        if regexp(name, fdrs{ifdr})
                            tokeep = false;
                            break
                        end
                    end
                    if tokeep
                        newpropfuncs{end + 1} = propfuncs{iprop};
                    end
                end
                propfuncs = newpropfuncs;

                
            end
            
        end
        
    end
end

