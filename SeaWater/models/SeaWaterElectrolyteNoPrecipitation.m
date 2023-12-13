classdef SeaWaterElectrolyteNoPrecipitation < ElectronicComponent

    properties

        K % Reaction rates for the reactions (depend on system)

        mainIonIndex

        con = PhysicalConstants(); % Physical constants 
        
        spdict     % dictionary with key={species names} and value={index in species concentration cell array}
        qpdict     % dictionary with key={quasi particle names} and value={index in quasi particle cell array}
        logspdict  % dictionary with key={main species names} and value={index in species log concentration cell array}
        solutedict % dictionary with key={solute names} and value={index in solute species cell array}
        soliddict  % dictionary with key={solid names} and value={index in solid species cell array}
        
        nsp     % number of species
        nlogsp  % number of main species (those given by logarithm)
        nqp     % number of quasi-particles
        nsolute % number of solute species
        nsolid  % number of solid species (those coming from precipitate, if any)
        
        indmainsp   % index mapping between main species and species index
        indsolutesp % index mapping between solute species index (solutedict) and global species index (spdict)        
        indspsolute % index mapping between global species index (spdict) and solute species index (solutedict). The index
                    % is zero is the given species is not a solute
        indsolidsp  % index mapping between solid species index (soliddict) and global species index (spdict)
        indspsolid  % index mapping between global species index (spdict) and solid species index (soliddict). The index
                    % is zero is the given species is not a solid        
        
        qpCompositionMatrix % composition matrix for the quasi particle:
                            % concentration(qp_i) = \sum_{j} qpCompositionMatrix_{i, j} concentration(sp_j)
        
        kappa % conductivity
        
        species % cell array for species 
                % Each cell is a struct with fields
                % - name (string)
                % - z (only for solutes, otherwise empty)
                % - lambda0 (only for solutes, otherwise empty)
                % - D : diffusion coefficients (if relevant, otherwise empty)
        
    end
    
    methods
        
        
        function model = SeaWaterElectrolyteNoPrecipitation(inputparams)
            
            model = model@ElectronicComponent(inputparams);
            
            fdnames = {'species'             , ...
                       'kappa'};
            model = dispatchParams(model, inputparams, fdnames);
            
            model = model.getAqueousMixtureComposition(inputparams);
            
            model = model.setupReactionRates(inputparams);
        
            model = model.setupMainIonIndex();

        end

        
        function model = setupReactionRates(model, inputparams)
            Kdata = inputparams.K;
            
            for ind = 1 : numel(Kdata)
                value  = Kdata(ind).value;
                str    = Kdata(ind).eval;
                K(ind) = eval(sprintf(str, value));
            end
            
            model.K = K;
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            nsp         = model.nsp;         % number of species
            nlogsp      = model.nlogsp;      % number of main species (those given by logarithm)
            nqp         = model.nqp;         % number of quasi-particles
            nsolute     = model.nsolute;     % number of solutes
            indmainsp   = model.indmainsp;   % index mapping between main species and species index
            indsolutesp = model.indsolutesp; % index mapping between solute species index (solutedict) and global species index (spdict)        
            indspsolute = model.indspsolute; % index mapping between global species index (spdict) and solute species index (solutedict)
        
            varnames = {};
            
            % electrolyte volumefraction
            varnames{end + 1} = 'volumeFraction'; 
            % species concentrations
            varnames{end + 1} = VarName({}, 'cs', nsp); 
            % species log-concentrations
            varnames{end + 1} = VarName({}, 'pcs', nlogsp);
            % species tranference numbers
            varnames{end + 1} = VarName({}, 'transNums', nsolute);
            % species tranference numbers (divided by concentratin to avoid singularity)
            varnames{end + 1} = VarName({}, 'transNumDcs', nsolute);
            % quasi particle total concentrations (with respect to total volume)
            varnames{end + 1} = VarName({}, 'qpepscs', nqp);
            % quasi particle total concentrations
            % varnames{end + 1} = VarName({}, 'qpcs', nqp);
            % quasi particle accumulation term
            varnames{end + 1} = VarName({}, 'qpAccums', nqp);
            % quasi particle diffusion flux
            varnames{end + 1} = VarName({}, 'qpDiffFluxes', nqp);
            % quasi particle migration flux
            varnames{end + 1} = VarName({}, 'qpMigFluxes', nqp);
            % quasi particle fluxes
            varnames{end + 1} = VarName({}, 'qpFluxes', nqp);
            % quasi particle source terms
            varnames{end + 1} = VarName({}, 'qpSrcTerms', nqp);
            % quasi particle mass conservations
            varnames{end + 1} = VarName({}, 'qpMassCons', nqp);
            % quasi particle mass conservations
            varnames{end + 1} = VarName({}, 'atomicMassCons', nqp);
            
            model = model.registerVarNames(varnames);
            
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateFromLogConcentration;
            for inlogsp = 1 : nlogsp
                inputnames = {VarName({}, 'pcs', nlogsp, inlogsp)};
                model = model.registerPropFunction({VarName({}, 'cs', nsp, indmainsp(inlogsp)), fn, inputnames});
            end
            
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateAtomicMassCons;
            for inqp = 1 : nqp
                inputnames = {VarName({}, 'qpepscs', nqp, inqp), VarName({}, 'cs', nsp), 'volumeFraction'};
                model = model.registerPropFunction({VarName({}, 'atomicMassCons', nqp, inqp), fn, inputnames});
            end
                        
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateDiffFluxes;
            inputnames = {VarName({}, 'cs', nsp), 'volumeFraction'};
            model = model.registerPropFunction({VarName({}, 'qpDiffFluxes', nqp), fn, inputnames});
            
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateTransferenceNumbers;
            inputnames = {VarName({}, 'cs', nsp, indsolutesp)};
            model = model.registerPropFunction({VarName({}, 'transNums', nsolute), fn, inputnames});           
            model = model.registerPropFunction({VarName({}, 'transNumDcs', nsolute), fn, inputnames});
            
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateMigFluxes;
            inputnames = {'j', VarName({}, 'transNums', nsolute)};
            model = model.registerPropFunction({VarName({}, 'qpMigFluxes', nqp), fn, inputnames});
            
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateFluxes;
            for inqp = 1 : nqp
                inputnames = {VarName({}, 'qpDiffFluxes', nqp, inqp), VarName({}, 'qpMigFluxes', nqp, inqp)};
                model = model.registerPropFunction({VarName({}, 'qpFluxes', nqp, inqp), fn, inputnames});
            end

            fn = @() SeaWaterElectrolyteNoPrecipitation.updateConductivity;
            inputnames = {'volumeFraction'};
            model = model.registerPropFunction({'conductivity', fn, inputnames});
            
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateCurrent;
            inputnames = {VarName({}, 'cs', nsp, indsolutesp), 'phi', 'conductivity', 'T', VarName({}, 'transNumDcs', nsolute)};
            model = model.registerPropFunction({'j', fn, inputnames});
            model = model.registerPropFunction({'jBcSource', fn, inputnames});
            
            
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateMassCons;
            for inqp = 1 : nqp
                inputnames = {VarName({}, 'qpAccums', nqp, inqp), ...
                              VarName({}, 'qpFluxes', nqp, inqp), ...
                              VarName({}, 'qpSrcTerms', nqp, inqp)};
                model = model.registerPropFunction({VarName({}, 'qpMassCons', nqp, inqp), fn, inputnames});
            end
            
            % update mass accumulations
            fn = @() SeaWaterElectrolyteNoPrecipitation.updateElectrolyteAccumTerms;
            fn = {fn, @(prop) PropFunction.accumFuncCallSetupFn(prop)};
            for iqp = 1 : nqp
                inputnames = {VarName({}, 'qpepscs', nqp, iqp)};
                model = model.registerPropFunction({VarName({}, 'qpAccums', nqp, iqp), fn, inputnames});
            end
            
            % Following commented functions are not used in assembly
            
            % fn = @() SeaWaterElectrolyteNoPrecipitation.computeLiquidQuasiParticleConcentration;
            % inputnames = {VarName({}, 'cs', nsp)};
            % model = model.registerPropFunction({VarName({}, 'qpcs', nqp), fn, inputnames});
            
        end
                
        function model = getAqueousMixtureComposition(model, inputparams)
            
            species        = arrayfun(@(sp) sp.name, inputparams.species, 'uniformoutput', false);
            quasiparticles = arrayfun(@(qp) qp.name, inputparams.quasiparticles, 'uniformoutput', false);
            solutes        = inputparams.solutes;
            logspecies     = inputparams.logspecies;
            solids         = inputparams.solids;
            
            nsp     = numel(species);
            nqp     = numel(quasiparticles);
            nlogsp  = numel(logspecies);            
            nsolute = numel(solutes);
            nsolid  = numel(solids);
            
            spdict     = containers.Map(species, 1 : nsp);
            qpdict     = containers.Map(quasiparticles, 1 : nqp);
            logspdict  = containers.Map(logspecies, 1 : nlogsp);
            solutedict = containers.Map(solutes, 1 : nsolute);
            soliddict  = containers.Map(solids, 1 : nsolid);
            
            keys = logspdict.keys;
            for ind = 1 : nlogsp
                indmainsp(logspdict(keys{ind})) = spdict(keys{ind});
            end
            
            keys = solutedict.keys;
            indspsolute = zeros(nsp, 1);
            for ind = 1 : nsolute
                indsolutesp(solutedict(keys{ind})) = spdict(keys{ind});
                indspsolute(spdict(keys{ind})) = solutedict(keys{ind});
            end

            keys = soliddict.keys;
            indspsolid = zeros(nsp, 1);
            for ind = 1 : nsolid
                indsolidsp(soliddict(keys{ind})) = spdict(keys{ind});
                indspsolid(spdict(keys{ind})) = soliddict(keys{ind});
            end

            
            qpmatrix = zeros(nqp, nsp);
            qp = inputparams.quasiparticles;
            for iqp = 1 : numel(qp)
                comp = qp(iqp).composition;
                for ic = 1 : numel(comp)
                    isp = spdict(comp(ic).name);
                    qpmatrix(iqp, isp) = comp(ic).coef;
                end
            end
            
            %% Return the values in the model
            
            model.qpCompositionMatrix = qpmatrix;

            model.spdict     = spdict;
            model.qpdict     = qpdict;
            model.logspdict  = logspdict;
            model.solutedict = solutedict;
            model.soliddict  = soliddict;
            
            model.indmainsp   = indmainsp;
            model.indsolutesp = indsolutesp;
            model.indspsolute = indspsolute;
            model.indsolidsp  = indsolidsp;
            model.indspsolid  = indspsolid;
            
            model.nsp     = nsp;
            model.nqp     = nqp;
            model.nlogsp  = nlogsp;
            model.nsolute = nsolute;
            model.nsolid  = nsolid;
            
        end

        function state = updateElectrolyteAccumTerms(model, state, state0, dt)

            vols = model.G.cells.volumes;
            nqp = model.nqp;

            for ind = 1 : nqp
                state.qpAccums{ind} =  1/dt*vols.*(state.qpepscs{ind} - state0.qpepscs{ind});
            end

        end
        
        function state = updateFromLogConcentration(model, state)
            
            nlogsp = model.nlogsp;
            indmainsp = model.indmainsp;
            
            pcs = state.pcs;
            
            for ind = 1 : nlogsp
                state.cs{indmainsp(ind)} = exp(pcs{ind});
            end
            
        end

        function state = updateAtomicMassCons(model, state)
            
            nqp         = model.nqp;
            C           = model.qpCompositionMatrix;
            indsolutesp = model.indsolutesp;
            
            qpepscs = state.qpepscs;
            cs      = state.cs;
            lvf     = state.volumeFraction; % liquid volume fraction
            
            for indqp = 1 : model.nqp
                atomicMassCons{indqp} = qpepscs{indqp};
                for isolute = 1 : model.nsolute
                    indsp = indsolutesp(isolute);
                    if C(indqp, indsp) ~= 0
                        atomicMassCons{indqp} = atomicMassCons{indqp} - C(indqp, indsp)*(lvf.*cs{indsp});
                    end
                end
            end
            
            state.atomicMassCons = atomicMassCons;
            
        end
        
        function state = updateDiffFluxes(model, state)
            species = model.species;
            nsolute     = model.nsolute;
            indsolutesp = model.indsolutesp;
            C = model.qpCompositionMatrix;
            
            cs = state.cs;
            vf = state.volumeFraction;
            effvf = vf.^1.5;
            
            for indqp = 1 : model.nqp
                qpDiffFluxes{indqp} = 0;
                for isolute = 1 : model.nsolute
                    indsp = indsolutesp(isolute);
                    if C(indqp, indsp) ~= 0
                        coef = (C(indqp, indsp)*species(indsp).D).*effvf;
                        qpDiffFluxes{indqp} = qpDiffFluxes{indqp} + assembleFlux(model, cs{indsp}, coef);
                    end
                end
            end
            
            state.qpDiffFluxes = qpDiffFluxes;
        end
        
        function state = updateTransferenceNumbers(model, state)
            
            nsolute     = model.nsolute;
            indsolutesp = model.indsolutesp;
            species     = model.species;
            
            cs = state.cs;
            
            total = 0;
            for ind = 1 : nsolute
                indsp = indsolutesp(ind);
                transNums{ind} = cs{indsp}.*abs(species(indsp).z).*species(indsp).lambda0;
                total = total + transNums{ind};
            end
            
            for ind = 1 : nsolute
                indsp = indsolutesp(ind);
                transNumDcs{ind} = abs(species(indsp).z).*species(indsp).lambda0./total;
                transNums{ind} = cs{indsp}.*transNumDcs{ind};
            end
            
            state.transNumDcs = transNumDcs;
            state.transNums   = transNums;
        end
        
        
        function state = updateMigFluxes(model, state)
            species = model.species;
            op      = model.operators;
            C       = model.qpCompositionMatrix;
            F       = model.con.F;
            
            j = state.j;
            transNums = state.transNums;
            
            for indqp = 1 : model.nqp
                migcoef_p = 0;
                migcoef_n = 0;
                for indsol = 1 : model.nsolute
                    indsp = model.indsolutesp(indsol);
                    if species(indsp).z > 0
                        migcoef_p = migcoef_p + C(indqp, indsp).*transNums{indsol}./(species(indsp).z.*F);
                    elseif species(indsp).z < 0
                        migcoef_n = migcoef_n + C(indqp, indsp).*transNums{indsol}./(species(indsp).z.*F);
                    end
                end
                qpMigFluxes{indqp} = (op.faceUpstr(j > 0, migcoef_p) + op.faceUpstr(j < 0, migcoef_n)).*j;
            end
            
            state.qpMigFluxes = qpMigFluxes;
            
        end

        function state = updateFluxes(model, state)
            
            nqp = model.nqp;
            
            qpMigFluxes = state.qpMigFluxes;
            qpDiffFluxes = state.qpDiffFluxes;
            
            for ind = 1 : nqp
                qpFluxes{ind} = qpMigFluxes{ind} + qpDiffFluxes{ind};
            end
            
            state.qpFluxes = qpFluxes;            
        end

        function state = updateConductivity(model, state)

            kappa = model.kappa;
            
            vf    = state.volumeFraction;
            
            state.conductivity = kappa.*vf.^1.5;
            
        end
        
        
        function state = updateCurrent(model, state)
            
            indsolutesp = model.indsolutesp;
            nsolute     = model.nsolute;
            R           = model.con.R;
            F           = model.con.F;
            species     = model.species;
            
            cs          = state.cs;
            phi         = state.phi;
            transNumDcs = state.transNumDcs;
            T           = state.T;
            kappa       = state.conductivity;
            % compute dmudc and jchems
            jchem = 0;
            dmu = R.*T;
            for ind = 1 : nsolute
                indsp = indsolutesp(ind);
                % dmudc = (R.*T)./cs{indsp}; to avoid singularities, we use transference number divided with
                % concentration (that is transNumDcs) and replace dmudc, which is equal to (R.*T)./concentration, with
                % dmu (dmu = RT, see above)
                if abs(species(indsp).z) > 0
                    jchemcoef = kappa.*transNumDcs{ind}.*dmu/(species(indsp).z*F);
                    jchem = jchem + assembleFlux(model, cs{indsp}, jchemcoef);
                end
            end
            
            j = assembleFlux(model, phi, kappa) + jchem;
            
            state.j = j;
            state.jBcSource = 0;
        end
        
        function state = updateMassCons(model, state)
            nqp = model.nqp;
            
            qpAccums   = state.qpAccums;
            qpFluxes   = state.qpFluxes;
            qpSrcTerms = state.qpSrcTerms;
            
            for ind = 1 : nqp
                qpMassCons{ind} = assembleConservationEquation(model, qpFluxes{ind}, 0, qpSrcTerms{ind}, qpAccums{ind});
            end
            
            state.qpMassCons = qpMassCons;
            
        end

                
        function state = computeLiquidQuasiParticleConcentration(model, state)
        
            nqp         = model.nqp;
            C           = model.qpCompositionMatrix;
            indsolutesp = model.indsolutesp;
            
            cs      = state.cs;
            
            for indqp = 1 : model.nqp
                qpcs{indqp} = 0;
                for isolute = 1 : model.nsolute
                    indsp = indsolutesp(isolute);
                    if C(indqp, indsp) ~= 0
                        qpcs{indqp} = qpcs{indqp} + C(indqp, indsp)*cs{indsp};
                    end
                end
            end
            
            state.qpcs = qpcs;
            
        end

        
        
    end
end
