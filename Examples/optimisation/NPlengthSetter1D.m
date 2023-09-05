classdef NPlengthSetter1D < LengthSetter1D

    properties

        NPratio % NP ratio
        alpha   % For a given NP ratio, alpha is equal the corresponding ratio between the thicknesses of the electrodes
        
    end
    
    methods
        
        function lengthsetter = NPlengthSetter1D(model, gridGenerator, NPratio)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            lengthsetter = lengthsetter@LengthSetter1D(gridGenerator, {ne, pe});
            lengthsetter.NPratio = NPratio;

            lengthsetter.alpha = lengthsetter.computeMultiplicationCoefficient(model);
            
        end
        
        function model = setLengths(lengthsetter, model, pelength)

            alpha = lengthsetter.alpha;
            
            v = double2ADI(zeros(2, 1), pelength);
            v(1) = alpha*pelength;
            v(2) = pelength;
            
            model = setLengths@LengthSetter1D(lengthsetter, model, v);
            
        end

        function v = getLengths(lengthsetter, model)

            v = lengthsetter.getAllLengths(model);
            v = v(3);

        end


        function alpha = computeMultiplicationCoefficient(lengthsetter, model)

            NPratio = lengthsetter.NPratio;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            
            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};

                theta100 = model.(elde).(am).(itf).theta100;
                theta0   = model.(elde).(am).(itf).theta0;
                cmax     = model.(elde).(am).(itf).cmax;
                vf       = unique(model.(elde).(am).volumeFraction);
                avf      = model.(elde).(am).activeMaterialFraction;
                F        = model.(elde).(am).constants.F;

                d.(elde) = abs(theta100 - theta0)*cmax*avf*vf*F/hour;

            end
            
            alpha = d.(pe)/d.(ne)*NPratio;
            
        end
    end
    
end
