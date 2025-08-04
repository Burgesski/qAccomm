function accommodationAnalysisFigures

    clear all;

    initialisationPlot = 0;
    animationCounter = 0; % Set to zero for no animation, or >1 for chron animation interval eg 10 to show every 10th iteration
    plotOneModelFlag = 1;
    shorelineInc = 10; % Cross section shorline markers might need a >1 time step interval to plot clearly
    lithoColourFlag = 1; % 0 = plot strata with water depth colour coding, 1 = plot output showing grainsize
    modelName = "clasticConstSLSSSandMud";
    params = parametersSiliciclasticConstSLSSSandMud(modelName, initialisationPlot);
    strata = qAccomm(modelName, params, animationCounter, plotOneModelFlag, lithoColourFlag);
    if ~animationCounter
        plotStrata(params, strata, params.totalChrons, modelName, animationCounter, lithoColourFlag, shorelineInc);
        plotModelAnalysis(params, strata, modelName)
    end

    initialisationPlot = 0;
    animationCounter = 0; % Set to zero for no animation, or >1 for chron animation interval eg 10 to show every 10th iteration
    plotOneModelFlag = 0;
    shorelineInc = 10; % Cross section shorline markers might need a >1 time step interval to plot clearly
    lithoColourFlag = 1; % 0 = plot strata with water depth colour coding, 1 = plot output showing grainsize
    modelName = "clasticVarSLSSSandMud";
    params = parametersSiliciclasticVarSLSSSandMud(modelName, initialisationPlot);
    strata = qAccomm(modelName, params, animationCounter, plotOneModelFlag, lithoColourFlag);
    if ~animationCounter
        plotStrata(params, strata, params.totalChrons, modelName, animationCounter, lithoColourFlag, shorelineInc);
        plotModelAnalysis(params, strata, modelName)
    end

    initialisationPlot = 0;
    animationCounter = 0; % Set to zero for no animation, or >1 for chron animation interval
    plotOneModelFlag = 0;
    shorelineInc = 1; % Set for consistency in function call, but not used in carbonate cross section plots
    lithoColourFlag = 0; % 0 = plot strata with water depth colour coding, 1 = plot output showing grainsize
    modelName = "carbonateVarSLSS";
    params = parametersCarbonateVarSeaSupplyRamp(modelName, initialisationPlot);
    strata = qAccomm(modelName, params, animationCounter, plotOneModelFlag, lithoColourFlag);
    if ~animationCounter
        plotStrata(params, strata, params.totalChrons, modelName, animationCounter, lithoColourFlag, shorelineInc);
        plotModelAnalysis(params, strata, modelName)
    end
end