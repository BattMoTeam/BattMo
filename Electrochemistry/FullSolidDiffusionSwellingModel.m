classdef FullSolidDiffusionSwellingModel < FullSolidDiffusionModel

    properties

    end

    methods

        function model = FullSolidDiffusionSwellingModel(paramobj)
            model = model@FullSolidDiffusionModel(paramobj);
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@FullSolidDiffusionModel(model);
                       
            fn = @FullSolidDiffusionSwellingModel.updateFlux;
            inputnames = {'c', 'D'};
            model = model.registerPropFunction({'flux', fn, inputnames});
            
            
            fn = @FullSolidDiffusionSwellingModel.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'radius', 'volumeFraction', 'Rvol'}});
            
            fn = @FullSolidDiffusionSwellingModel.updateMassAccum;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccum', fn, {'c'}});
            
            fn = @FullSolidDiffusionSwellingModel.assembleSolidDiffusionEquation;
            model = model.registerPropFunction({'solidDiffusionEq', fn, {'c', 'cSurface', 'massSource', 'D'}});

            
        end
        

        function state = updateMassSource(model, state)
            %%Modification of mass source
           
            op  = model.operators;
            rp  = state.radius;
            rp0 = model.rp;
            vf  = state.volumeFraction;
            amf = model.activeMaterialFraction;
            
            Rvol = state.Rvol;

            Rvol = op.mapFromBc*Rvol;


            L = length(Rvol.val);


            %just to define a generic AD structure for massSource
            A = model.rp;
            state.massSource = - Rvol*((4*pi*A^3)/(3*amf*A));

            % Same expression as in non-swelling material but taking into
            % account the variation of radius and volume fraction as well
            % as a normalisation term at the end as we are in fixed coordinates)
            for i = 1:L
                R = Rvol.val(i);
                n = fix(i/L);
                radius = rp.val(n+1);
                volumeFraction = vf.val(n+1);
                state.massSource.val(i) = - R*((4*pi*radius^3)/(3*amf*volumeFraction))*(rp0/radius)^3;
            end


            
        end

        function state = updateDiffusionCoefficient(model, state)

            if model.useDFunc

                computeD = model.computeDFunc;
                cmax     = model.cmax;
                theta0   = model.theta0;
                theta100 = model.theta100;
                
                c = state.c;

                cmin = theta0*cmax;
                cmax = theta100*cmax;

                soc = (c - cmin)./(cmax - cmin);
                
                D = computeD(soc);

                state.D = D;
                
            else
                
                state = updateDiffusionCoefficient@SolidDiffusionModel(model, state);
                
            end

        end

        function state = updateMassAccum(model, state, state0, dt)

            op = model.operators;      
            c = state.c;
            c0 = state0.c;
            
            massAccum = 1/dt*op.vols.*(c - c0);
            state.massAccum = massAccum;
            
        end
        
        function state = updateMassConservation(model, state)
           
            op = model.operators;
            
            flux       = state.flux;
            massSource = state.massSource;
            massAccum  = state.massAccum;

            
            state.massCons = massAccum + op.div(flux) - massSource;

        end
        
        function state = updateFlux(model, state)
            
            useDFunc = model.useDFunc;

            op = model.operators;
            c = state.c;
            D = state.D;

            if useDFunc
                state.flux = op.flux(D, c);
            else
                D = op.mapToParticle*D;
                state.flux = op.flux(D, c);
            end
            
            
        end
    
        function state = assembleSolidDiffusionEquation(model, state)
            
        %% TODO : change name of this function
            
            op       = model.operators;
            useDFunc = model.useDFunc;
            
            c     = state.c;
            D     = state.D;
            cSurf = state.cSurface;
            src   = state.massSource;

            if ~useDFunc
                % TODO : make this implementation better
                % Here, we first dispatch D on all the particle cells and, then, retain only the value at the boundary.
                D = op.mapToParticle*D;
            end
            
            D = op.mapToBc*D;
            
            eq = D.*op.Tbc.*(op.mapToBc*c - cSurf) + op.mapToBc*src;
            
            state.solidDiffusionEq = eq;
            
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
    
