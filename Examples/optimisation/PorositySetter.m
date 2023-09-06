classdef PorositySetter

    properties
        
        nxs        % number of cells for each component 
        compnames   % names of the components
        compinds   % index of the components
        elytemaps
        refvalues
        
    end
    
    methods
        
        function porosetter = PorositySetter(model, components)

            ne    = 'NegativeElectrode';
            sep   = 'Separator';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            am    = 'ActiveMaterial';
            
            compnames = {ne, sep, pe};

            compinds = ismember({ne, sep, pe}, components);

            elyteinds = zeros(model.(elyte).G.parentGrid.getNumberOfCells(), 1);
            elyteinds(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.getNumberOfCells())';
            
            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                switch compname
                  case ne
                    compmodel = model.(ne).(am);
                  case sep
                    compmodel = model.(elyte).(sep);
                  case pe
                    compmodel = model.(pe).(am);
                  otherwise
                    error('compname not recognized');
                end
                nxs.(compname)  = compmodel.G.getNumberOfCells();
                poro = unique(compmodel.porosity);
                if (numel(poro) > 1)
                    warning('The initial model does not have homogeneous porosity');
                end
                refvalues(icomp) = poro;
                elytemaps.(compname) = elyteinds(compmodel.G.mappings.cellmap);
                
            end
            
            porosetter.nxs       = nxs;
            porosetter.compnames = compnames;
            porosetter.compinds  = compinds;
            porosetter.refvalues = refvalues;
            porosetter.elytemaps = elytemaps;
            
        end
        
        function model = setPorosities(porosetter, model, v)

            nxs       = porosetter.nxs;
            compnames = porosetter.compnames;
            compinds  = porosetter.compinds;
            refvalues = porosetter.refvalues;
            
            poros = refvalues;
            poros = subsasgnAD(poros, compinds, v);

            for icomp = 1 : numel(compnames)
                compname = compnames{icomp};
                model = porosetter.setComponentPorosity(model, poros(icomp), compname);
            end

        end

        function model = setComponentPorosity(porosetter, model, poro, compname)
        % We update the porosities and the dependent values.

            nxs       = porosetter.nxs;
            elytemaps = porosetter.elytemaps;
            
            v = poro.*ones(nxs.(compname), 1);

            ne    = 'NegativeElectrode';
            sep   = 'Separator';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            am    = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            switch compname
              case {ne, pe}
                model.(compname).(am).porosity             = v;
                model.(compname).(am).volumeFraction       = 1 - v;
                switch model.(compname).(am).diffusionModelType
                  case {'simple', 'interParticleOnly'}
                    % nothing to do
                  case 'full'
                    model.(compname).(am).(sd).volumeFraction = 1 - v;
                  otherwise
                    error('Unknown diffusionModelType %s', diffusionModelType);
                end
              case sep
                model.(elyte).(sep).porosity       = v;
                model.(elyte).(sep).volumeFraction = 1 - v;
              otherwise
                error('compname not recognized');
            end
            model.(elyte).volumeFraction = subsasgnAD(model.(elyte).volumeFraction, elytemaps.(compname), v);

        end
            
        
        function v = getAllPorosities(porosetter, model)

            compnames = porosetter.compnames;
            
            ne    = 'NegativeElectrode';
            sep   = 'Separator';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            am    = 'ActiveMaterial';
            
            v = nan(numel(compnames), 1);

            for icomp = 1 : numel(compnames)
                
                compname = compnames{icomp};
                switch compname
                  case {ne, pe}
                    poro = model.(compname).(am).porosity;
                  case sep
                    poro = model.(elyte).(sep).volumeFraction;
                  otherwise
                    error('compname not recognized');
                end
                
                v = subsasgnAD(v, icomp, poro(1));
                
            end
            
        end


        function v = getPorosities(porosetter, model)

            compinds = porosetter.compinds;

            v = porosetter.getAllPorosities(model);
            v = v(compinds);
            
        end
    end
    
end
