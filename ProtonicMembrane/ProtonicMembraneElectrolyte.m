classdef ProtonicMembraneElectrolyte < BaseModel
    
    properties
        
        % Temperature
        T
        % Structure with physical constants
        constants
        % Equilibrium H2 potential
        EH2_0
        % Equilibrium O2 potential
        EO2_0
        % Equibrium p-type conductivity 
        sigmaP_0        
        % Equibrium n-type conductivity 
        sigmaN_0
        % proton conductivity
        sigmaHp

        d_micro = 30 % thickness in micron

        Y
        dH_hyd  % kJ/mol
        dS_hyd  % J / mol / K
        Ea_prot % Activation energy proton diffusion
        pH2O_in
        pH2O_neg
        SU

        t_p_O2  % tr.nr holes in 1bar humidified oxygen
        
    end
    
    methods
        
        function model = ProtonicMembraneElectrolyte(paramobj)

            model = model@BaseModel();
            
            fdnames = {'G'       , ...
                       'T'       , ...
                       'EH2_0'   , ...
                       'EO2_0'   , ...
                       'sigmaN_0', ...
                       'Y'       , ...
                       'dH_hyd'  , ...
                       'dS_hyd'  , ...
                       'Ea_prot' , ...
                       'pH2O_in' , ...
                       'pH2O_neg', ...
                       'SU'      , ...
                       't_p_O2' };
            
            model = dispatchParams(model, paramobj, fdnames);

            model.operators = localSetupOperators(model.G);

            model.constants = PhysicalConstants();

            con = model.constants;

            % Compute pressures
            
            pH2O = model.pH2O_in*(1 - model.SU);

            pH2 = model.pH2O_in;

            pO2 = (1 - model.SU)*model.pH2O_in/2;
            pO2 = 1;

            % Compute reaction constants
            
            Y_mol   = model.Y/((4.22e-8)^3*con.Na);
            D0_prot = 0.021*38/model.T;       % pre-exp proton diffusion
            D_prot  = D0_prot*exp(-(model.Ea_prot*1000)/con.R/model.T);
            
            K_hyd   = exp((model.dS_hyd/con.R) - model.dH_hyd*1000/(con.R*model.T));
            K_H     = K_hyd*pH2O;
            K_H_neg = K_hyd*model.pH2O_neg;
            
            Y = model.Y;
            OH_pos  = ((3*K_H - sqrt(K_H*(9*K_H - 6*K_H*Y + K_H*Y^2 + 24*Y - 4*Y^2)))/(K_H - 4));  % Per formula unit
            OH_neg = ((3*K_H_neg - sqrt(K_H_neg*(9*K_H_neg - 6*K_H_neg*Y + K_H_neg*Y^2 + 24*Y - 4*Y^2)))/(K_H_neg - 4));

            % Computate sigmaHp

            sigma_prot_pos = (con.F*OH_pos*D_prot);
            sigma_prot_neg = (con.F*OH_pos*D_prot);

            sigmaHp = (sigma_prot_pos + sigma_prot_neg)/2;

            % Computate sigmaP_0

            t_p_O2  = 0.5; % tr.nr holes in 1bar humidified oxygen
            K_H_p   = K_hyd*0.027;
            p_ref   = ((3*K_H_p - sqrt(K_H_p*(9*K_H_p - 6*K_H_p*Y + K_H_p*Y^2 + 24*Y - 4*Y^2)))/(K_H_p - 4)); 
            sigmaP_0 = (t_p_O2/(1 - t_p_O2))*(con.F*p_ref*D_prot); 

            
            % Assign the values to the model
            
            model.sigmaHp  = sigmaHp;
            model.sigmaP_0 = sigmaP_0;
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % Electrostatic potential phi
            varnames{end + 1} = 'phi';
            % Electromotive potential pi (see Jacobsen and Mogensen paper) also called Electrochemical potential
            varnames{end + 1} = 'pi';
            % Electronic Chemical potential E
            varnames{end + 1} = 'E';
            % Electronic conductivity
            varnames{end + 1} = 'sigmaEl';
            % Hp conductivity
            varnames{end + 1} = 'sigmaHp';
            % H+ source term
            varnames{end + 1} = 'sourceHp';            
            % Electronic source term
            varnames{end + 1} = 'sourceEl';
            % H+ flux
            varnames{end + 1} = 'jHp';            
            % Electronic flux
            varnames{end + 1} = 'jEl';
            % H+ mass conservation (we measure the mass in Coulomb, hence "massConsHp")
            varnames{end + 1} = 'massConsHp';
            % Charge conservation
            varnames{end + 1} = 'chargeConsEl';            
            % Charge conservation
            varnames{end + 1} = 'alpha';
            
            model = model.registerVarNames(varnames);

            fn = @ProtonicMembraneElectrolyte.updateE;
            inputnames = {'phi', 'pi'};
            model = model.registerPropFunction({'E', fn, inputnames});
            
            fn = @ProtonicMembraneElectrolyte.updateSigmaHp;
            inputnames = {};
            model = model.registerPropFunction({'sigmaHp', fn, inputnames});

            fn = @ProtonicMembraneElectrolyte.updateHpFlux;
            inputnames = {'phi'};
            model = model.registerPropFunction({'jHp', fn, inputnames});

            fn = @ProtonicMembraneElectrolyte.updateElConductivity;
            inputnames = {'E', 'alpha'};
            model = model.registerPropFunction({'sigmaEl', fn, inputnames});
                        
            fn = @ProtonicMembraneElectrolyte.updateElFlux;
            inputnames = {'sigmaEl', 'pi'};
            model = model.registerPropFunction({'jEl', fn, inputnames});

            fn = @ProtonicMembraneElectrolyte.updateMassConsHp;
            inputnames = {'sourceHp', 'jHp'};
            model = model.registerPropFunction({'massConsHp', fn, inputnames});
        
            fn = @ProtonicMembraneElectrolyte.updateChargeConsEl;
            inputnames = {'sourceEl', 'jEl'};
            model = model.registerPropFunction({'chargeConsEl', fn, inputnames});
        
        end

        function state = updateE(model, state)

            state.E = state.phi - state.pi;
            
        end

        function state = updateSigmaHp(model, state)

            nc = model.G.cells.num;
           
            state.sigmaHp = model.sigmaHp*ones(nc, 1);
            
        end
        
        function state = updateHpFlux(model, state)

            op = model.operators;
            
            sigmaHp = state.sigmaHp;
            phi     = state.phi;
            
            state.jHp = assembleFlux(model, phi, sigmaHp);

        end

        function state = updateElConductivity(model, state)

            F = model.constants.F;
            R = model.constants.R;
            T = model.T;
            
            sigmaP_0 = model.sigmaP_0;
            sigmaN_0 = model.sigmaN_0;
            EH2_0    = model.EH2_0;
            EO2_0    = model.EO2_0;
            
            E     = state.E;
            alpha = state.alpha;
            
            f = F/(R*T);
            sigmaEl = sigmaP_0*exp(f*(E - EO2_0)) + sigmaN_0*exp(-f*(E - EH2_0));
            sigmaEl = (1 - alpha)*(sigmaP_0 + sigmaN_0) + alpha*sigmaEl;
            
            state.sigmaEl = sigmaEl;
            
        end
        
        function state = updateElFlux(model, state)

            sigmaEl = state.sigmaEl;
            pi      = state.pi;
            
            state.jEl = assembleFlux(model, pi, sigmaEl);
            
        end

        
        function state = updateChargeConsEl(model, state)
            
            op = model.operators;
            
            sourceEl = state.sourceEl;
            jEl      = state.jEl;
            
            state.chargeConsEl =  op.Div(jEl) - sourceEl;
            
        end
        
        function state = updateMassConsHp(model, state)
            
            op = model.operators;
            
            sourceHp = state.sourceHp;
            jHp      = state.jHp;
            
            state.massConsHp =  op.Div(jHp) - sourceHp;
            
        end
        
    end
    
end
