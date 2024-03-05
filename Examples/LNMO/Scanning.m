classdef Scanning < handle

    properties (SetAccess = immutable)
        
        jsonstruct
        
    end

    properties

        model
        specificEnergyHour % in Wh/hour
        specificEnergy % 
        thickness 
        SiContent

    end
    
    methods

        function sc = Scanning()

            jsonstruct      = parseBattmoJson('Examples/LNMO/SiGrLnmo.json');
            jsonstruct_cccv = parseBattmoJson('Examples/JsonDataFiles/cccv_control.json');

            jsonstruct = removeJsonStructFields(jsonstruct, ...
                                                {'Control', 'DRate'}        , ...
                                                {'Control', 'controlPolicy'}, ...
                                                {'Control', 'dEdtLimit'}    , ...
                                                {'Control', 'dIdtLimit'}    , ...
                                                {'Control', 'rampupTime'}   , ...
                                                {'Control', 'useCVswitch'});

            jsonstruct = mergeJsonStructs({jsonstruct, jsonstruct_cccv});


            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            sep = 'Separator';
            co  = 'Coating';
            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            jsonstruct.(ne).(co).(am1).(itf).guestStoichiometry100 = 0.99;
            jsonstruct.(ne).(co).(am1).(itf).guestStoichiometry0   = 0.01;
            jsonstruct.(ne).(co).(am2).(itf).guestStoichiometry100 = 0.99;
            jsonstruct.(ne).(co).(am2).(itf).guestStoichiometry0   = 0.01;

            jsonstruct.(pe).(co).(am).(itf).guestStoichiometry0   = 0.99;
            jsonstruct.(pe).(co).(am).(itf).guestStoichiometry100 = 0.01;

            % We do not need fine discretization

            jsonstruct.(ne).(co).N = 1;
            jsonstruct.(pe).(co).N = 1;
            jsonstruct.(sep).N = 1;

            sc.jsonstruct = jsonstruct;
            
        end

        function computeSpecificEnergy(sc, thickness, SiContent)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            sep = 'Separator';
            co  = 'Coating';
            bd  = 'Binder';
            ad  = 'ConductingAdditive';
            am1 = 'ActiveMaterial1';
            am2 = 'ActiveMaterial2';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            jsonstruct = sc.jsonstruct;
            
            jsonstruct.(ne).(co).thickness = thickness;
            jsonstruct.(ne).(co).(am2).massFraction = SiContent;
            jsonstruct.(ne).(co).(am1).massFraction = 1 - (jsonstruct.(ne).(co).(am2).massFraction + ...
                                                           jsonstruct.(ne).(co).(bd).massFraction + ...
                                                           jsonstruct.(ne).(co).(ad).massFraction);
            
            model = setupModelFromJson(jsonstruct);

            [capacity, capacities] = computeCellCapacity(model);

            npratio = capacities.(ne)/capacities.(pe);

            jsonstruct.(pe).(co).thickness = jsonstruct.(pe).(co).thickness*(npratio/1.1);

            model = setupModelFromJson(jsonstruct);
            
            energy = computeCellEnergy(model);
            mass   = computeCellMass(model);

            specificEnergy = energy/mass;

            sc.specificEnergy     = specificEnergy;
            sc.specificEnergyHour = specificEnergy/hour;
            sc.thickness          = thickness;
            sc.SiContent          = SiContent;
            sc.model              = model;
        end
        
    end

    
end

