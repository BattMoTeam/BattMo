setupoutput = false;
if setupoutput
    ind = cellfun(@(state) ~isempty(state), states);
    states = states(ind);
    for istate = 1 : numel(states)
        states{istate} = model.addVariables(states{istate});
    end
end

% set(0, 'defaultlinelinewidth', 1);
% set(0, 'defaultaxesfontsize', 15);


close all
figure
hold on
for istate = 1 : numel(states)
    % val = states{istate}.(oer).(ptl).liqeps;
    % val = states{istate}.(oer).(ctl).etaElyte;
    % val = states{istate}.(her).(ptl).phasePressures{1};
    % val = states{istate}.(her).(ptl).compGasSources{2};
    % val = states{istate}.(her).(ptl).eSource;
    % val = states{istate}.(her).(ptl).H2OliquidSource;
    % val = states{istate}.(her).(ptl).OHsource;
    % val = states{istate}.(her).(ptl).concentrations{1};
    % val = states{istate}.(her).(ptl).liqrhoeps;
    % val = states{istate}.(her).(ptl).liquidMassFlux;
    % val = states{istate}.(her).(ptl).liqrho;
    % val = states{istate}.(her).(ptl).phasePressures{1};
    % val = states{istate}.(her).(ptl).viscosities{1};
    % val = states{istate}.(her).(ptl).phaseFluxes{1};
    % val = states{istate}.(her).(ctl).etaInmr;
    % val = states{istate}.(her).(ctl).inmrReactionRate;
    % val = states{istate}.(her).(ctl).inmrOHsource;
    % val = states{istate}.(her).(exl).OHexchangeRate;
    % val = states{istate}.(her).(exl).H2OexchangeRate;
    % val = states{istate}.(her).(exl).H2OaElyte;
    val = states{istate}.(her).(ptl).liqeps;
    % val = states{istate}.(her).(ptl).conductivity;
    % val = states{istate}.(inm).conductivity;
    % val = states{istate}.(inm).eSource;
    % val = states{istate}.(inm).OHsource;
    % val = states{istate}.(inm).phi;
    % val = states{istate}.(inm).H2Oceps;
    % val = states{istate}.(inm).H2Oc;
    % val = states{istate}.(inm).H2OdiffFlux;
    % val = states{istate}.(inm).H2OmigFlux;
    % val = states{istate}.(inm).H2OSource;
    h = plot(val);
    if istate == numel(states)
        set(h, 'linewidth', 3);
    end
    % val = states{istate}.(her).(exl).H2OaInmr;
    % val = states{istate}.(her).(ptl).OHsource;
    % h = plot(val);
    % if istate == numel(states)
        % set(h, 'linewidth', 2);
    % end
end

return

%%

%%
coupcells = model.couplingTerms{1}.couplingcells;

figure

clear esource
inmresource = states{1}.(inm).eSource;
esource(coupcells(:, 1)) = inmresource(coupcells(:, 2));

plot(esource)

figure

clear j
inmrj = states{1}.(inm).j;
j(coupcells(1 : end - 1, 1)) = inmrj(coupcells(1 : end - 1, 2));

plot(j)

figure

clear chargecons
inmrchargecons = states{1}.(inm).chargeCons;
chargecons(coupcells(1 : end, 1)) = inmrchargecons(coupcells(1 : end, 2));

plot(chargecons)
