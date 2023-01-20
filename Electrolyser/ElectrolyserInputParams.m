classdef ElectrolyserInputParams < InputParams
    
    properties

        G % parent grid (the component grids are subgrid of that one)

        IonomerMembrane
        HydrogenEvolutionElectrode
        OxygenEvolutionElectrode        
                
        couplingTerms
        
        controlI % given value for galvanistic control
    end
    
    methods
        
        function paramobj = ElectrolyserInputParams(jsonstruct)

            paramobj = paramobj@InputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.IonomerMembrane = IonomerMembraneInputParams(pick('IonomerMembrane'));
            paramobj.HydrogenEvolutionElectrode = EvolutionElectrodeInputParams(pick('HydrogenEvolutionElectrode'));
            paramobj.OxygenEvolutionElectrode = EvolutionElectrodeInputParams(pick('OxygenEvolutionElectrode'));

            paramobj = paramobj.validateInputParams();
        end

        function paramobj = validateInputParams(paramobj)

            assert(strcmp(paramobj.HydrogenEvolutionElectrode.porousTransportLayerType, 'Hydrogen'), 'Expected porous transport layer is hydrogen');
            assert(strcmp(paramobj.OxygenEvolutionElectrode.porousTransportLayerType, 'Oxygen'), 'Expected porous transport layer is oxygen');

            paramobj = validateInputParams@InputParams(paramobj);
            
        end
        
    end
    
    
end
