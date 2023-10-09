function [energy, dischargeFunction] = computeCellEnergy(model, varargin)
%% Compute the cell theoritcal maximum energy.
%% The discharge function is the voltage as a function of the state of charge for an infinitely slow discharge.
    
    opt = struct('capacities' , [], ...
                 'temperature', 298);
    opt = merge_options(opt, varargin{:});

    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    am  = 'ActiveMaterial';
    co  = 'Coating';
    itf = 'Interface';
    sd  = 'SolidDiffusion';

    if isempty(opt.CRate)
        
        if isempty(opt.capacities)
            [~, capacities] = computeCellCapacity(model);
        else
            capacities = opt.capacities
        end

        T = opt.temperature;

        capacity = min(capacities.(ne), capacities.(pe));


        th100 = 'guestStoichiometry100';
        th0   = 'guestStoichiometry0';
        sc    = 'saturationConcentration';
        
        % setup concentration at start and end of discharge for each electrode
        itfmodel = model.(ne).(co).(am).(itf);
        c0s{1} = itfmodel.(th100)*itfmodel.(sc);
        cTs{1} = itfmodel.(th0)*itfmodel.(sc);
        itfmodel = model.(pe).(co).(am).(itf);
        c0s{2} = itfmodel.(th0)*itfmodel.(sc);
        cTs{2} = itfmodel.(th100)*itfmodel.(sc);

        N = 1000;

        eldes = {ne, pe};
        
        for ielde = 1 : numel(eldes)

            elde = eldes{ielde};
            
            smax = capacity./capacities.(elde);
            
            s = smax.*linspace(0, 1, N + 1)';

            c    = (1 - s).*c0s{ielde} + s.*cTs{ielde};
            cmax = model.(elde).(co).(am).(itf).(sc);
            
            f = model.(elde).(co).(am).(itf).computeOCPFunc(c(1 : end - 1), T, cmax);
            
            % function handler
            fs{ielde} = @(s) model.(elde).(co).(am).(itf).computeOCPFunc((1 - s).*c0s{ielde} + s.*cTs{ielde}, T, cmax);
            
            energies{ielde} = capacities.(elde)*smax/N*sum(f);
            
        end
        
        energy = (energies{2} - energies{1});

        dischargeFunction = @(s) (fs{2}(s) - fs{1}(s));

    else

        CRate = opt.CRate;
        [energy, dischargeFunction] = computeCellEnergyGivenCrate(model, CRate, extra{:});
        
    end
    
end
