function [energy, output] = computeCellEnergy(model, varargin)
%% Compute the cell energy.
%
% If no CRate is given, it corresponds to the maximum theoritical cell energy, at infinitly small CRate. Otherwise,
% simulate the discharge and gives the corresponding energy.
%
% The output structure provides more detailed information with the fields
% - energy
% - dischargeFunction % function handler giving the voltage as a function of the state of charge
% - E                 % Voltage output (raw computation data in case CRate is given)
% - I                 % Current output (raw computation data in case CRate is given)
% - time              % time output (raw computation data in case CRate is given)
    
    opt = struct('capacities' , [], ...
                 'temperature', 298, ...
                 'CRate', []);
    [opt, extra] = merge_options(opt, varargin{:});

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

        N = 1000;

        eldes = {ne, pe};
        
        for ielde = 1 : numel(eldes)

            elde = eldes{ielde};
            
            smax = capacity./capacities.(elde);
            
            itfmodel = model.(elde).(co).(am).(itf);
            c0 = itfmodel.(th100)*itfmodel.(sc);
            cT = itfmodel.(th0)*itfmodel.(sc);            
            
            s = smax.*linspace(0, 1, N + 1)';

            c    = (1 - s).*c0 + s.*cT;
            cmax = model.(elde).(co).(am).(itf).(sc);
            
            f = model.(elde).(co).(am).(itf).computeOCPFunc(c(1 : end - 1), T, cmax);
            
            % function handler
            fs{ielde} = @(s) model.(elde).(co).(am).(itf).computeOCPFunc((1 - s).*c0 + s.*cT, T, cmax);
            
            energies{ielde} = capacities.(elde)*smax/N*sum(f);
            
        end
        
        energy = (energies{2} - energies{1});

        dischargeFunction = @(s) (fs{2}(s) - fs{1}(s));

        output = struct('energy'           , energy, ...
                        'dischargeFunction', dischargeFunction);
        
    else

        CRate = opt.CRate;
        output = computeCellEnergyGivenCrate(model, CRate, extra{:});
        
    end

    energy = output.energy;
    
end
