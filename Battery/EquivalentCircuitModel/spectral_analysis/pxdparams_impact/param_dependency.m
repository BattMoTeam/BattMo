function param_dependency()
    mrstModule add ad-core mrst-gui mpfa agmg linearsolvers
    
    % 1. Chargement des données (Conservez votre méthode de chargement ici)
    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
    jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});
    
    c_ne = 29.866*mol/litre; 
    c_pe = 17.038*mol/litre; 
    
    % 2. Création de l'interface graphique (Fenêtre plus haute pour accueillir les curseurs)
    fig = uifigure('Name', 'Impedance Explorer', 'Position', [100, 50, 600, 650]);
    
    % Graphique
    ax = uiaxes(fig, 'Position', [50, 350, 500, 250]);
    title(ax, 'Nyquist Diagram');
    xlabel(ax, 'Re(Z)');
    ylabel(ax, '-Im(Z)');
    grid(ax, 'on');
    
    % Label de statut
    lbl = uilabel(fig, 'Position', [50, 310, 500, 22], 'Text', 'Ready, use slider to start.');
    lbl.FontWeight = 'bold';

    % --- CURSEUR 1 : Diffusion Solide (Anode) ---
    % Modifie le transport du Lithium dans les particules de graphite.
    uilabel(fig, 'Position', [50, 260, 500, 40], 'WordWrap', 'on', ...
        'Text', '1. refDiffusionCoefficient * 10^x');

    sld1 = uislider(fig, 'Position', [50, 240, 500, 3], ...
        'Limits', [-10, 10], ...           % Définit le minimum et le maximum
        'MajorTicks', -10:1:10, ...        % Ajoute une graduation tous les 1
        'Value', 0);
    
    % --- CURSEUR 2 : Taux de réaction (Anode) ---
    % Modifie la cinétique de transfert de charge.
    uilabel(fig, 'Position', [50, 170, 500, 40], 'WordWrap', 'on', ...
        'Text', '2. Anode reactionRateConstant * 10^x');
    sld2 = uislider(fig, 'Position', [50, 150, 500, 3],...
        'Limits', [-10, 10], ...           % Définit le minimum et le maximum
        'MajorTicks', -10:1:10, ...        % Ajoute une graduation tous les 1
        'Value', 0);
    
    % --- CURSEUR 3 : Taux de réaction (Cathode) ---
    % Modifie la cinétique de transfert de charge.
    uilabel(fig, 'Position', [50, 80, 500, 40], 'WordWrap', 'on', ...
        'Text', '3. Cathode reactionRateConstant * 10^x');
    sld3 = uislider(fig, 'Position', [50, 60, 500, 3], ...
        'Limits', [-10, 10], ...           % Définit le minimum et le maximum
        'MajorTicks', -10:1:10, ...        % Ajoute une graduation tous les 1
        'Value', 0);

    % 4. Connexion aux événements
    % On crée une fonction anonyme qui passe TOUS les curseurs
    callback_fcn = @(src, event) updatePlot(ax, lbl, jsonstruct, c_ne, c_pe, sld1, sld2, sld3);
    
    sld1.ValueChangedFcn = callback_fcn;
    sld2.ValueChangedFcn = callback_fcn;
    sld3.ValueChangedFcn = callback_fcn;
    
    % Optionnel : Tracer l'état initial
    % updatePlot(ax, lbl, jsonstruct, c_ne, c_pe, sld1, sld2, sld3);
end

% --- Fonction de mise à jour ---
function updatePlot(ax, lbl, jsonstruct, c_ne, c_pe, sld1, sld2, sld3)
    
    % Figer l'affichage
    lbl.Text = 'Processing... ';
    lbl.FontColor = '#D95319';
    drawnow; 
    
    % 1. Récupération des valeurs de base depuis votre JSON
    base_diff_anode = 1.3135e-15;
    base_rate_anode = 5.031e-11;
    base_rate_cathode = 2.33e-11;
    
    % 2. Application des multiplicateurs (10^slider_val)
    val_diff_anode = base_diff_anode * (10^(sld1.Value));
    val_rate_anode = base_rate_anode * (10^(sld2.Value));
    val_rate_cathode = base_rate_cathode * (10^(sld3.Value));
    
    % 3. Injection dans la structure JSON AVANT de construire le modèle
    jsonstruct.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.referenceDiffusionCoefficient = val_diff_anode;
    jsonstruct.NegativeElectrode.Coating.ActiveMaterial.Interface.reactionRateConstant = val_rate_anode;
    jsonstruct.PositiveElectrode.Coating.ActiveMaterial.Interface.reactionRateConstant = val_rate_cathode;
    
    try
        % 4. Exécution de BattMo
        [model, inputparams, ~] = setupModelFromJson(jsonstruct);
        initstate = initStateChen2020(model, c_ne, c_pe);
        
        impsolv = ImpedanceSolver(inputparams, 'initstate', initstate, 'computeSteadyState', false);
    
        frequences = logspace(-2, 3, 30); 
        Z = impsolv.computeImpedance(frequences);
        
        % 5. Tracé
        plot(ax, real(Z), -imag(Z), '-o', 'LineWidth', 2, 'Color', '#0072BD');
        
        % Retour à la normale
        lbl.Text = 'Done';
        lbl.FontColor = '#77AC30';
        
    catch ME
        lbl.Text = 'Erreur lors du calcul (vérifiez la console).';
        lbl.FontColor = '#A2142F';
        disp('--- ERREUR DANS LE CALCUL ---');
        disp(ME.message);
        for k = 1:length(ME.stack)
            disp(['Ligne ', num2str(ME.stack(k).line), ' dans ', ME.stack(k).name]);
        end
    end
end