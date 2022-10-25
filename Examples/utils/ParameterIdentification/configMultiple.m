function config = configMultiple(jsonExp)

    configBV     = configButlerVolmer(jsonExp);
    configBrugge = configBruggeman(jsonExp);
    configArea   = configVolumetricSurfaceArea(jsonExp);

    config = [configBV;
              configBrugge;
              configArea];

end
