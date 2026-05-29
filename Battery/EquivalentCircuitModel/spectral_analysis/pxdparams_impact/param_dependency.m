function param_dependency()
    mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
    
    jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                                   jsonstruct_geometry});
    
    
    
    includeDoubleLayer = false;
    
    if includeDoubleLayer
    
        jsonstruct.(ne).(co).(am).(itf).useDoubleLayerCapacity = true;
        jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = 0.2;
    
    end
    
    
    [~, inputparams, ~] = setupModelFromJson(jsonstruct);
    
    c_ne = 29.866*mol/litre; % initial concentration at negative electrode
    c_pe = 17.038*mol/litre; % initial concentration at positive electrode
    
    % 2. Création de l'interface graphique
    fig = uifigure('Name', 'Impedance explorer', 'Position', [100, 100, 600, 450]);
    ax = uiaxes(fig, 'Position', [50, 150, 500, 250]);
    title(ax, 'Diagramme de Nyquist (Chen 2020)');
    xlabel(ax, 'Z'' (\Omega)');
    ylabel(ax, '-Z" (\Omega)');
    grid(ax, 'on');

    % 3. Création du curseur (-1 à +1 ordre de grandeur)
    sld = uislider(fig, ...
        'Position', [100, 80, 400, 3], ...
        'Limits', [-1, 1], ... 
        'Value', 0);
    
    lbl = uilabel(fig, 'Position', [100, 100, 400, 22], 'Text', 'Prêt. Bougez le curseur pour lancer un calcul.');

    % 4. Connexion au relâchement de la souris (ValueChangedFcn)
    % On passe la structure inputparams complète à notre fonction
    sld.ValueChangedFcn = @(src, event) updatePlot(ax, lbl, inputparams, event.Value);
    
    % (Optionnel) Lancer un premier calcul d'initialisation
    % updatePlot(ax, lbl, inputparams, 0);
end

% --- Fonction de mise à jour lourde ---
function updatePlot(ax, lbl, inputparams, slider_val)
    
    % 1. Avertir l'utilisateur et figer l'affichage
    lbl.Text = 'Calcul en cours... Veuillez patienter (cela peut prendre quelques secondes).';
    drawnow; % CRUCIAL : Force MATLAB à afficher le texte avant de bloquer
    
    % 2. Calculer la nouvelle valeur et modifier la structure
    % On modifie ici le SolidDiffusion.referenceDiffusionCoefficient
    valeur_initiale = 3.3e-14;
    nouvelle_valeur = valeur_initiale * (10^slider_val);
    
    % On navigue dans l'arborescence exacte de votre fichier JSON
    inputparams.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.referenceDiffusionCoefficient = nouvelle_valeur;
    
    % 3. APPEL DE VOS FONCTIONS
    try
        [model, ~, ~] = setupModelFromJson(jsonstruct);
        initstate = initStateChen2020(model, c_ne, c_pe);
        impsolv = ImpedanceSolver(inputparams, 'initstate', initstate, 'computeSteadyState', false);
    
        % Définition des fréquences et calcul de l'impédance
        frequences = logspace(-2, 3, 30); 
        Z = impsolv.computeImpedance(frequences);
        
        % 4. Tracer le vrai résultat
        % On utilise -imag(Z) car on trace traditionnellement l'opposé de la partie imaginaire
        plot(ax, real(Z), -imag(Z), '-o', 'LineWidth', 2, 'Color', '#D95319');
        
        % Restaurer le texte
        lbl.Text = sprintf('Calcul terminé. Diffusion = %.2e (10^{%.2f}x)', nouvelle_valeur, slider_val);
        
    catch ME
        % Si le solveur plante (par exemple si la nouvelle valeur rend le système instable)
        lbl.Text = 'Erreur lors du calcul (paramètre hors limites ?)';
        disp(ME.message); % Affiche l'erreur dans la console MATLAB
    end
end