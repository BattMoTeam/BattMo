classdef SwellingMaterial < ActiveMaterial
    
    properties

    end
    
    methods
        
        function model = SwellingMaterial(paramobj)
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
            model = model@ActiveMaterial(paramobj);
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ActiveMaterial(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {{sd, 'cAverage'} ,...
                        'radius'   , ...
                        'porosity' ,...
                        'volumeFraction'};
            model = model.registerVarNames(varnames);
      

            fn = @SwellingMaterial.updatePorosity;
            model = model.registerPropFunction({'porosity', fn, {}});

            fn = @SwellingMaterial.updateEffectiveElectricalConductivity;
            model = model.registerPropFunction({'EffectiveElectricalConductivity', fn, {'volumeFraction'} });


            if model.use_particle_diffusion
                
                fn = @SwellingMaterial.updateRadius;
                model = model.registerPropFunction({'radius', fn, {{sd, 'cAverage'}}});
                model = model.registerPropFunction({{sd, 'radius'}, fn, {{sd, 'cAverage'}} });
                
                fn = @SwellingMaterial.updateVolumeFraction;
                model = model.registerPropFunction({{'volumeFraction'}, fn, {'porosity'}});
                model = model.registerPropFunction({{sd, 'volumeFraction'}, fn, {'porosity'}});
                model = model.registerPropFunction({{itf, 'volumeFraction'}, fn, {'porosity'}});


                fn = @SwellingMaterial.updateVolumetricSurfaceArea;
                model = model.registerPropFunction({{itf, 'volumetricSurfaceArea'}, fn, {'radius', {itf, 'volumeFraction'}}});
                model = model.registerPropFunction({{sd, 'volumetricSurfaceArea'}, fn, {{sd, 'radius'},{sd, 'volumeFraction'}}});

                
            else

                fn = @SwellingMaterial.updateRadius;
                model = model.registerPropFunction({'radius', fn, {{sd, 'cAverage'}}});
                
                fn = @ActiveMaterial.updateVolumeFraction;
                model = model.registerPropFunction({'volumeFraction', fn, {'porosity'}});
                model = model.registerPropFunction({{itf, 'volumeFraction'}, fn, {'porosity'}});

                fn = @SwellingMaterial.updateVolumetricSurfaceArea;
                model = model.registerPropFunction({{itf, 'volumetricSurfaceArea'}, fn, {{sd,'radius'}, {itf, 'volumeFraction'}}});

                
            end


            fn = @SwellingMaterial.SolidDiffusionModel.updateAverageConcentration;
            % TODO check definition of update function for cAverage in FullSolidDiffusionModel (updateAverageConcentration)
            % If it can be used, then Xavier fix that
            model = model.registerPropFunction({{sd, 'cAverage'}, fn, {{sd, 'c'}}});

            
                      
                   
        end

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);

            state = updateControl(model, state, drivingForces);


            %new
            state                = model.updatePorosity(state, state0, dt);
            state                = model.updateVolumeFraction(state);
            state                = model.updateEffectiveElectricalConductivity(state)    ;     
           
           
            state                = model.updateStandalonejBcSource(state);
            state                = model.updateCurrent(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassAccum(state.SolidDiffusion, state0.SolidDiffusion, dt);
            state                = model.dispatchTemperature(state);
            state.SolidDiffusion = model.SolidDiffusion.updateDiffusionCoefficient(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateFlux(state.SolidDiffusion);
            state                = model.updateConcentrations(state);
            
            %new
            state                = model.updateAverageConcentration(state);
            state                = model.updateRadius(state);
            state                = model.updateVolumetricSurfaceArea(state);
            state                = model.updateOperators(state.SolidDiffusion);

            state                = model.updatePhi(state);
            state.Interface      = model.Interface.updateReactionRateCoefficient(state.Interface);
            state.Interface      = model.Interface.updateOCP(state.Interface);
            state.Interface      = model.Interface.updateEta(state.Interface);
            state.Interface      = model.Interface.updateReactionRate(state.Interface);
            state                = model.updateRvol(state);
            state                = model.updateCurrentSource(state);
            state                = model.updateChargeConservation(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassSource(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.assembleSolidDiffusionEquation(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateMassConservation(state.SolidDiffusion);
            
            %% Setup equations and add some scaling
            n     = model.(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F     = model.(sd).constants.F;
            vol   = model.operators.pv;
            rp    = state.radius;
            vsf   = state.Interface.volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
            
            eqs = {};
            eqs{end + 1} = state.chargeCons;
            eqs{end + 1} = scalingcoef*state.(sd).massCons;
            eqs{end + 1} = scalingcoef*state.(sd).solidDiffusionEq;
            
            names = {'chargeCons', ...
                     'massCons', ...
                     'solidDiffusionEq'};
            
            types = {'cell', 'cell', 'cell'};

            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end


        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@ActiveMaterial(model, state, problem, dx, drivingForces);
            
        end





        function state = updateConcentrations(model, state)

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            if model.use_particle_diffusion
            if strcmp(model.diffusionModelType, 'simple')
                state.(sd).cAverage = state.c;
            end
            
            state.(itf).cElectrodeSurface = state.(sd).cSurface;
            else
                state.(itf).cElectrodeSurface = state.c;
            end
            
        end

        function state = updateMassFlux(model, state)
        % Used when diffusionModelType == 'simple'

            D = model.EffectiveDiffusionCoefficient;
            
            c = state.c;

            massflux = assembleFlux(model, c, D);
            
            state.massFlux = massflux;

        end
            
        function state = assembleAccumTerm(model, state, state0, dt)
        % Used when diffusionModelType == 'simple'
            
            vols   = model.G.cells.volumes;
            vf     = state.volumeFraction;
            amFrac = model.activeMaterialFraction;

            c  = state.c;
            c0 = state0.c;

            state.massAccum = vols.*vf.*amFrac.*(c - c0)/dt;
            
        end


        
        function state = updateMassConservation(model, state)
        % Used when diffusionModelType == 'simple' or no particle diffusion
            
            flux = state.massFlux;
            source = state.massSource;
            accum = state.massAccum;

            cons = assembleConservationEquation(model, flux, 0, source, accum);
            
            state.massCons = cons;
            
        end
        
        function state = updateStandalonejBcSource(model, state)
            
            state.jBcSource = state.controlCurrentSource;

        end

        function state = updateCurrentSource(model, state)
            
            F    = model.Interface.constants.F;
            vols = model.G.cells.volumes;
            n    = model.Interface.n;

            Rvol = state.Rvol;
            
            state.eSource = - vols.*Rvol*n*F; % C/s
            
        end
        
        function state = updatePhi(model, state)
            state.Interface.phiElectrode = state.phi;
        end         
        
        function state = dispatchTemperature(model, state)
            state.Interface.T = state.T;
            state.SolidDiffusion.T = state.T;
        end

        function state = updatejBcSource(model, state)
            state.jBcSource = state.jCoupling + state.jExternal;
        end
        
        function state = updatejFaceBc(model, state)
            state.jFaceBc = state.jFaceCoupling + state.jFaceExternal;
        end
        
        function state = updatejExternal(model, state)
            state.jExternal = 0;
            state.jFaceExternal = 0;
        end

        function state = updatejCoupling(model, state)
            state.jCoupling = 0;
            state.jFaceCoupling = 0;
        end


        function state = updateAverageConcentrationSD(model, state)
            %replaces the function updateAverageConcentration of the
            %FullSolidDiffusion class

            % shortcut
            sd  = 'SolidDiffusion';

            op = state.(sd).operators;
            vols = op.vols;
            map = op.mapToParticle;
            
            c = state.(sd).c;

            m    = map'*(c.*vols); % total amount [mol] in the cell particles
            vols = map'*(vols);    % volume 

            state.(sd).cAverage = m./vols;
            
        end
        

        function state = updateAverageConcentration(model, state)

            % shortcut
            sd  = 'SolidDiffusion';

            vf       = state.volumeFraction;
            am_frac  = model.activeMaterialFraction;
            vols     = model.G.cells.volumes;
            
            c = state.(sd).cAverage;

            vols = am_frac*vf.*vols;

            cAverage = sum(c.*vols)/sum(vols);

            state.cAverage = cAverage;
            
        end
        
        
        function state = updateSOC(model, state)

            % shortcut
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            vf       = state.volumeFraction;
            am_frac  = model.activeMaterialFraction;
            vols     = model.G.cells.volumes;
            cmax     = model.(itf).cmax;
            theta100 = model.(itf).theta100;
            theta0   = model.(itf).theta0;
            
            c = state.(sd).cAverage;

            theta = c/cmax;
            m     = (1 ./ (theta100 - theta0));
            b     = -m .* theta0;
            SOC   = theta*m + b;
            vol   = am_frac*vf.*vols;
            
            SOC = sum(SOC.*vol)/sum(vol);

            state.SOC = SOC;
            
        end
   

%Functions specific to swelling materials
        function state = updateRadius(model, state)

            c = state.(sd).cAverage;
            
            radius_0 = model.SolidDiffusion.rp;
            densitySi = model.density;
            MolarMassSi = 28.0855 * 1E-3;
            molarVolumeSi = MolarMassSi/densitySi;
            molarVolumeLi = 8.8 * 1E-6;
            cmaxLi = model.Interface.cmax;
            
            radius = radius_0 * (1 + (3.75*molarVolumeLi*c)/(cmaxLi*molarVolumeSi))^(1/3);
            state.radius = radius;

            if model.use_particle_diffusion
                state.SolidDiffusion.radius = radius;
            end
         end

            
        function state = updatePorosity(model, state, state0, dt)
            porosity = state0.porosity;

            a = state0.volumetricSurfaceArea;       
            R = state0.Interface.R;

            molarVolumeLithiated = model.updateMolarVolumeLithiated(state);
            densitySi = model.Interface.density;
            molarMassSi = 28.0855 * 1E-3;
            molarVolumeSi = molarMassSi/densitySi;

            state.porosity = porosity + dt*a*R*(molarVolumeLithiated - 3.75*molarVolumeSi);
        end


        function state = updateVolumeFraction(model, state)
            porosity = state.porosity;
            vf = 1 - porosity;

            state.volumeFraction = vf;
            state.Interface.volumeFraction = vf;

             if model.use_particle_diffusion
                state.SolidDiffusion.volumeFraction = vf;
            end
        end
           
       
        function state = updateVolumetricSurfaceArea(model, state)
            vf = state.Interface.volumeFraction;
            amf = model.activeMaterialFraction;
            radius = state.SolidDiffusion.radius;

            vsa = 3*vf*amf/radius;
            state.Interface.volumetricSurfaceArea = vsa;

            if model.use_particle_diffusion
                state.SolidDiffusion.volumetricSurfaceArea = vsa;
            end
        end


        function state = updateEffectiveElectricalConductivity(model, state)
            state.EffectiveElectricalConductivity; 
            vf = 1 - state.porosity;
            brugg = model.BruggemanCoefficient;
                   
            % setup effective electrical conductivity using Bruggeman approximation
            state.EffectiveElectricalConductivity = model.electricalConductivity.*vf.^brugg;
        end
      
                
        function molarVolumeLitihated = updateMolarVolumeLithiated(model, state)

            sd  = 'SolidDiffusion';

            c = state.(sd).cAverage;

            densitySi = model.Interface.density;
            molarMassSi = 28.0855 * 1E-3;
            molarVolumeSi = molarMassSi/densitySi;

            molarVolumeLi = 8.8 * 1E-6;
            cmaxLi = model.Interface.cmax;

            molarVolumeLitihated = 4/15*(molarVolumeSi + 3.75*(c/cmaxLi)*molarVolumeLi);
        end
        
        function state = updateOperators(model, state)
            np = model.SolidDiffusion.np;
            
            N  = model.SolidDiffusion.N;
            rp = state.radius;
            
            celltbl.cells = (1 : np)';
            celltbl = IndexArray(celltbl);

            % Solid particle cells
            Scelltbl.Scells = (1 : N)';
            Scelltbl = IndexArray(Scelltbl);

            cellScelltbl = crossIndexArray(celltbl, Scelltbl, {}, 'optpureproduct', true);
            cellScelltbl = sortIndexArray(cellScelltbl, {'cells', 'Scells'});

            endScelltbl.Scells = N;
            endScelltbl = IndexArray(endScelltbl);
            endcellScelltbl = crossIndexArray(cellScelltbl, endScelltbl, {'Scells'});

            G = cartGrid(N, rp); 
            r = G.nodes.coords;

            G.cells.volumes   = 4/3*pi*(r(2 : end).^3 - r(1 : (end - 1)).^3);
            G.cells.centroids = (r(2 : end) + r(1 : (end - 1)))./2;

            G.faces.centroids = r;
            G.faces.areas     = 4*pi*r.^2;
            G.faces.normals   = G.faces.areas;

            rock.perm = ones(N, 1);
            rock.poro = ones(N, 1);

            tbls = setupSimpleTables(G);
            cellfacetbl = tbls.cellfacetbl;

            hT = computeTrans(G, rock); % hT is in cellfacetbl

            cells = cellfacetbl.get('cells');
            faces = cellfacetbl.get('faces');
            sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1; % sgn is in cellfacetbl
                
            % We change name of cellfacetbl
            ScellSfacetbl = cellfacetbl;
            ScellSfacetbl = replacefield(ScellSfacetbl, {{'cells', 'Scells'}, {'faces', 'Sfaces'}});
            
            % Here, we use that we know *apriori* the indexing in G.cells.faces (the last index corresponds to outermost cell-face)
            Tbc = hT(end); % half-transmissibility for of the boundary face
            Tbc = repmat(Tbc, np, 1);
            
            Sfacetbl.Sfaces = (2 : N)'; % index of the internal faces (correspond to image of C')
            Sfacetbl = IndexArray(Sfacetbl);
            cellSfacetbl = crossIndexArray(celltbl, Sfacetbl, {}, 'optpureproduct', true);

            % we consider only the internal faces
            allScellSfacetbl = ScellSfacetbl;
            ScellSfacetbl = crossIndexArray(allScellSfacetbl, Sfacetbl, {'Sfaces'});

            cellScellSfacetbl = crossIndexArray(celltbl, ScellSfacetbl, {}, 'optpureproduct', true);

            map = TensorMap();
            map.fromTbl = allScellSfacetbl;
            map.toTbl = cellScellSfacetbl;
            map.mergefds = {'Scells', 'Sfaces'};
            map = map.setup();            

            hT = map.eval(hT);
            sgn = map.eval(sgn);

            %% setup of divergence operator

            prod = TensorProd();
            prod.tbl1 = cellScellSfacetbl;
            prod.tbl2 = cellSfacetbl;
            prod.tbl3 = cellScelltbl;
            prod.mergefds = {'cells'};
            prod.reducefds = {'Sfaces'};
            prod = prod.setup();

            divMat = SparseTensor();
            divMat = divMat.setFromTensorProd(sgn, prod);
            divMat = divMat.getMatrix();

            div = @(u) (divMat*u);
            
            gradMat = -divMat';

            prod = TensorProd();
            prod.tbl1 = cellScellSfacetbl;
            prod.tbl2 = cellScelltbl;
            prod.tbl3 = cellSfacetbl;
            prod.mergefds = {'cells'};
            prod.reducefds = {'Scells'};
            prod = prod.setup();

            invHtMat = SparseTensor();
            invHtMat = invHtMat.setFromTensorProd(1./hT, prod);
            invHtMat = invHtMat.getMatrix();

            flux = @(D, c) -(1./(invHtMat*(1./D))).*(gradMat*c);

            %% External flux map (from the boundary conditions)

            map = TensorMap();
            map.fromTbl = endcellScelltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'cells', 'Scells'};
            map = map.setup();

            f = map.eval(ones(endcellScelltbl.num, 1));

            prod = TensorProd();
            prod.tbl1 = cellScelltbl;
            prod.tbl2 = celltbl;
            prod.tbl3 = cellScelltbl;
            prod.mergefds = {'cells'};

            mapFromBc = SparseTensor();
            mapFromBc = mapFromBc.setFromTensorProd(f, prod);
            mapFromBc = mapFromBc.getMatrix();

            mapToBc = mapFromBc';

            %% map from cell (celltbl) to cell-particle (cellScelltbl)
            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'cells'};
            map = map.setup();

            mapToParticle = SparseTensor();
            mapToParticle = mapToParticle.setFromTensorMap(map);
            mapToParticle = mapToParticle.getMatrix();
            
            vols = G.cells.volumes;

            map = TensorMap();
            map.fromTbl = Scelltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'Scells'};
            map = map.setup();

            vols = map.eval(vols);

            state.SolidDiffusion.operators = struct('div'          , div          , ...
                               'flux'         , flux         , ...
                               'mapFromBc'    , mapFromBc    , ...
                               'mapToParticle', mapToParticle, ...
                               'mapToBc'      , mapToBc      , ...
                               'Tbc'          , Tbc          , ...
                               'vols'         , vols);
            
        end


        function state = updateFlux(model, state)
            
            useDFunc = model.SolidDiffusion.useDFunc;
            op = state.SolidDiffusion.operators;
            
            c = state.SolidDiffusion.c;
            D = state.SolidDiffusion.D;

            if useDFunc
                state.SolidDiffusion.flux = op.flux(D, c);
            else
                D = op.mapToParticle*D;
                state.SolidDiffusion.flux = op.flux(D, c);
            end          
        end


        function state = updateMassSource(model, state)
            if strcmp(model.diffusionModelType, 'simple')
            % used when diffusionModelType == simple
            
            vols = model.G.cells.volumes;
            Rvol = state.Rvol;
            state.massSource = - Rvol.*vols;

            else
            % used when diffusionModelType == full
            op  = state.SolidDiffusion.operators;
            rp  = state.SolidDiffusion.radius;
            vf  = model.SolidDiffusion.volumeFraction;
            amf = model.SolidDiffusion.activeMaterialFraction;
            Rvol = state.SolidDiffusion.Rvol;

            Rvol = op.mapFromBc*Rvol;
            
            state.SolidDiffusion.massSource = - Rvol*((4*pi*rp^3)/(3*amf*vf));
            end
            
        end
       

        function state = updateMassAccum(model, state, state0, dt)

            op = state.SolidDiffusion.operators;
            
            c = state.c;
            c0 = state0.c;
            
            state.SolidDiffusion.massAccum = 1/dt*op.vols.*(c - c0);
            
        end
        






        
    end
    
end

%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
