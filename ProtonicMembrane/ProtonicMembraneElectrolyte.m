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
        
    end
    
    methods
        
        function model = ProtonicMembraneElectrolyte(paramobj)

            model = model@BaseModel();
            
            fdnames = {'G'};
            model = dispatchParams(model, paramobj, fdnames);

            model.operators = localSetupOperators(model.G);

            model.constants = PhysicalConstants();
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % Electrostatic potential phi
            varnames{end + 1} = 'phi';
            % Electromotive potential pi (see Jacobsen and Mogensen paper) or Electrochemical potential
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
            % H+ mass conservation (we measure the mass in Coulomb, hence "chargeConsHp")
            varnames{end + 1} = 'chargeConsHp';
            % Charge conservation
            varnames{end + 1} = 'chargeConsEl';            

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
            inputnames = {'E'};
            model = model.registerPropFunction({'sigmaEl', fn, inputnames});
                        
            fn = @ProtonicMembraneElectrolyte.updateElFlux;
            inputnames = {'sigmaEl', 'pi'};
            model = model.registerPropFunction({'jEl', fn, inputnames});

            fn = @ProtonicMembraneElectrolyte.updateChargeConsHp;
            inputnames = {'sourceHp', 'jHp'};
            model = model.registerPropFunction({'chargeConsHp', fn, inputnames});
        
            fn = @ProtonicMembraneElectrolyte.updateChargeConsEl;
            inputnames = {'sourceEl', 'jEl'};
            model = model.registerPropFunction({'chargeConsEl', fn, inputnames});
        
        end
        
        
        function state = updateHpFlux(model, state)

            sigmaHp = model.sigmaHp;
            op = model.operators;
            
            phi = state.phi;
            
            state.jHp = -sigmaHp*op.Grad(phi);
            
        end
        
        function state = updateElConductivity(model, state)
            
            F = model.constants.F;
            R = model.constants.R;
            T = model.T;
            
            sigmaP_0 = model.sigmaP_0;
            sigmaN_0 = model.sigmaN_0;
            EH2_0    = model.EH2_0;
            EO2_0    = model.EO2_0;
            
            E = state.E;

            f = F/(R*T);
            
            sigmaEl = sigmaP_0*exp(-f*(E - EO2_0)) + sigmaN_0*exp(-f*(E - EH2_0));
            
            state.sigmaEl = sigmaEl;
            
        end
        
        function state = updateElFlux(model, state)
            
            sigmaEl = state.sigmaEl;
            phi     = state.phi;
            E       = state.E;
            
            state.jEl = assembleFlux(model, phi + E, sigmaEl);
            
        end
        
        function state = updateChargeConsEl(model, state)
            
            op = model.operators;
            
            sourceEl = state.sourceEl;
            jEl      = state.jEl;
            
            state.chargeConsEl =  op.div(jEl) - sourceEl;
            
        end
        
        function state = updateChargeConsHp(model, state)
            
            op = model.operators;
            
            sourceEl = state.sourceHp;
            jHp      = state.jHp;
            
            state.chargeConsHp =  op.div(jHp) - sourceHp;
            
        end
        
    end
    
end
