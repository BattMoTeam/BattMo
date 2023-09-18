function [energy, dischargeFunction] = computeCellEnergy(model, varargin)
%% Compute the cell theoritcal maximum energy.
%% The discharge function is the voltage as a function of the state of charge for an infinitely slow discharge.
    
    opt = struct('capacities' , [], ...
                 'temperature', 298, ...
                 'CRate', []);
    [opt, extra] = merge_options(opt, varargin{:});

    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    am  = 'ActiveMaterial';
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

        % setup concentration at start and end of discharge for each electrode
        itfmodel = model.(ne).(am).(itf);
        c0s{1} = itfmodel.theta100*itfmodel.cmax;
        cTs{1} = itfmodel.theta0*itfmodel.cmax;
        itfmodel = model.(pe).(am).(itf);
        c0s{2} = itfmodel.theta0*itfmodel.cmax;
        cTs{2} = itfmodel.theta100*itfmodel.cmax;

        N = 1000;

        eldes = {ne, pe};
        
        for ielde = 1 : numel(eldes)

            elde = eldes{ielde};
            
            smax = capacity./capacities.(elde);
            
            s = smax.*linspace(0, 1, N + 1)';

            c    = (1 - s).*c0s{ielde} + s.*cTs{ielde};
            cmax = model.(elde).(am).(itf).cmax;
            
            f = model.(elde).(am).(itf).computeOCPFunc(c(1 : end - 1), T, cmax);
            
            % function handler
            fs{ielde} = @(s) model.(elde).(am).(itf).computeOCPFunc((1 - s).*c0s{ielde} + s.*cTs{ielde}, T, cmax);
            
            energies{ielde} = capacities.(elde)*smax/N*sum(f);
            
        end
        
        energy = (energies{2} - energies{1});

        dischargeFunction = @(s) (fs{2}(s) - fs{1}(s));

    else

        CRate = opt.CRate;
        [energy, dischargeFunction] = computeCellEnergyGivenCrate(model, CRate, extra{:});
        
    end
    
end
