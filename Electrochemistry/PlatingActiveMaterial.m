classdef PlatingActiveMaterial < ActiveMaterial

    methods

        function model = PlatingActiveMaterial(inputparams)
            model = model@ActiveMaterial(inputparams);
            if model.useLithiumPlating
                model.LithiumPlating = LithiumPlatingLatz();
            end
        end

        function model = registerVarAndPropfuncNames(model)
            model = registerVarAndPropfuncNames@ActiveMaterial(model);
        end

        function state = updateRvol(model, state)
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            lp  = 'LithiumPlating';

            vsa = model.(itf).volumetricSurfaceArea;
            R = state.(itf).R;

            if model.useLithiumPlating
                theta = state.(lp).surfaceCoverage;
                Rchem = state.(lp).chemicalFlux;
                Rvol = vsa .* (R .* (1 - theta) + Rchem .* theta);
            else
                Rvol = vsa .* R;
            end

            state.(sd).Rvol = Rvol;
        end

    end

end