classdef OxygenElectrode < OpenElectrode
    

    properties
        sp.O2.MW
    end
    
    
    
    methods
        
        function model = OxygenElectrode(paramobj)

        % paramobj is instance of ElectronicComponentInputParams
            model = model@OpenElectrode(paramobj);
            
            model.compInd.O2 = 5;
            model.gasInd.O2 = 2;
            
        end

        function state = updateGasPressure(model, state)
            MWH2O = model.sp.H2O.MW;
            MWO2 = model.sp.O2.MW;
            
            mO2  = state.O2rhoeps;
            mH2O = state.H2Ogasrhoeps;
            vfg  = state.volumeFractions{model.phaseInd.gas};
            

            state.pressures{model.phaseInd.gas} = (mH20./MWH2O + mO2./MWO2).*R.*T./vfg;
        end
        
        
        function state = updateGasViscosity(model, state)
            T = state.T;

            state.viscosities(model.phaseInd.gas) = (0.1971 + T * (0.0803 - 3.99e-5 * T)) * 1e-6;
        end
        
        
    end


end

