classdef ProtonicMembraneAnode < ProtonicMembraneElectrode
    
    properties
        
        % coefficient in Buttler-Volmer
        beta
        
        % charge-transfer current density
        i_0

        % Limiting current densities
        ila % Anode
        ilc % Cathode

        n
        
        R_ct_0
        Ea_ct      
        SU         
        O2_conc_feed
        steam_ratio
        Ptot
        
        pO2
        pH2O
        pRef
        
    end
    
    methods
        
        function model = ProtonicMembraneAnode(paramobj)

            model = model@ProtonicMembraneElectrode(paramobj);

            fdnames = {'beta'        , ...
                       'ila'         , ...
                       'ilc'         , ...
                       'R_ct_0'      , ...
                       'Ea_ct'       , ...
                       'n'           , ...
                       'SU'          , ...
                       'O2_conc_feed', ...
                       'steam_ratio' , ...
                       'Ptot'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            model.constants = PhysicalConstants();

            model.pRef = 1*barsa;
            
            c = model.constants;

            T = model.T;
            
            % compute pressures to setup Eocv
            
            pO2_in  = model.O2_conc_feed*(1 - model.steam_ratio)*model.Ptot;
            pH2O_in = model.steam_ratio*model.Ptot;

            pO2  = pO2_in + model.SU*pH2O_in/2;
            pH2O = pH2O_in*(1 - model.SU);

            pRef = model.pRef;
            
            % Compute charge-transfer current density
            
            R_ct = model.R_ct_0.*exp(model.Ea_ct./(c.R.*T)).*(pO2/pRef).^(-0.2); 
            f = c.F/(c.R*T);
            i_0 = 1./(R_ct.*f.*model.n); % charge-transfer current density

            % Setup the gas components

            if strcmp(model.gasSupplyType, 'notCoupled')
                model.gasSupplyType = 'coupled'
                fprintf('We switch to gasSupplyType to coupled in Anode model\n');
            end
            
            gasInd.H2O = 1;
            gasInd.O2  = 2;
            nGas = 2;
            
            % Assign the computed values
            model.i_0    = i_0;
            model.nGas   = nGas;
            model.gasInd = gasInd;
            model.pO2    = pO2;
            model.pH2O   = pH2O;
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@ProtonicMembraneElectrode(model);

            varnames = {};
            % charge Transfer Current Density
            varnames{end + 1} = 'chargeTransferCurrentDensity';

            model = model.registerVarNames(varnames);
            
            fn = @ProtonicMembraneAnode.updateJHp;
            inputnames = {'eta', 'chargeTransferCurrentDensity'};
            model = model.registerPropFunction({'jHp', fn, inputnames});

            fn = @ProtonicMembraneAnode.updateI0;
            inputnames = {VarName({}, 'pressures', model.nGas)};
            model = model.registerPropFunction({'chargeTransferCurrentDensity', fn, inputnames});

            fn = @ProtonicMembraneAnode.updatePressures;
            inputnames = {};
            model = model.registerPropFunction({VarName({}, 'pressures', model.nGas), fn, inputnames});

            fn = @ProtonicMembraneAnode.updateOcp;
            inputnames = {VarName({}, 'pressures', model.nGas)};
            model = model.registerPropFunction({'Eocv', fn, inputnames});
            
            % fn = @ProtonicMembraneAnode.updateEta2;
            % inputnames = {'j', 'alpha'};
            % model = model.registerPropFunction({'eta', fn, inputnames});

            % fn = @ProtonicMembraneAnode.updatePhi2;
            % inputnames = {'eta', 'pi', 'Eocv'};
            % model = model.registerPropFunction({'phi', fn, inputnames});

        end

        function state = updateOcp(model, state)

            c      = model.constants;
            gasInd = model.gasInd;
            T      = model.T;
            pRef   = model.pRef;

            pH2O = state.pressures{gasInd.H2O};
            pO2  = state.pressures{gasInd.O2};
            
            state.Eocv = model.E_0 - 0.00024516.*T - c.R.*T./(2*c.F)*log((pH2O/pRef)./((pO2/pRef).^(1/2)));
            
        end
        
        function state = updatePressures(model, state)

            gasInd = model.gasInd;
            N = model.N;

            onevec = ones(N, 1);
            state.pressures{gasInd.H2O} = model.pH2O*onevec;
            state.pressures{gasInd.O2}  = model.pO2*onevec;
            
        end

        function state = updateI0(model, state)
            
            c      = model.constants;
            gasInd = model.gasInd;
            T      = model.T;
            pRef   = model.pRef;
            
            pO2 = state.pressures{gasInd.O2};
            
            R_ct = model.R_ct_0.*exp(model.Ea_ct./(c.R.*T)).*(pO2/pRef).^(-0.2); 
            f = c.F/(c.R*T);
            i0 = 1./(R_ct.*f.*model.n);
            
            state.chargeTransferCurrentDensity = i0;
            
        end
        
        function state = updateJHp(model, state)
            
            con  = model.constants;

            beta = model.beta;
            ila  = model.ila;
            ilc  = model.ilc;
            
            eta = state.eta;
            i0  = state.chargeTransferCurrentDensity;
            
            feta = con.F*model.n/(con.R*model.T).*eta;
            
            jHp = -i0.*(exp(-beta*feta) - exp((1 - beta)*feta))./(1 + (i0./ilc).*exp(-beta*feta) - (i0./ila).*exp((1 - beta)*feta));

            % jHp = i0*(exp(feta/2) - exp(-feta/2))/2;
            
            % R = 5;
            % jHp = 1/R*eta;
            
            state.jHp = jHp;
            
        end

        % function state = updateEta2(model, state)

        %     R = 0.05;
            
        %     j = state.j;
        %     eta = R*j; 
            
        %     state.eta = eta;
            
        % end

        % function state = updatePhi2(model, state)

        %     state.phi = state.pi - state.eta - state.Eocv;
            
        % end

        
    end
    
end

