function params = parametersCarbonateVarSeaSupplyRamp(modelName, plotFlag)

    rng(42); % Standard value 42

    params = struct; % Create params structure array, primarily to make sure it is empty
    params.siliciclasticModel = 0;
    params.carbonateModel = 1;
    params.totalChrons = 1000; % total number of chron surfaces calculated in model time steps
    params.chronInterval = 0.001; % Duration of each chron My
    params.AATimeStep = 0.00001; % Duration of each time step in the accommodation analysis (AA) calculation
    params.mcIterations = 100; % Number of realizations calculated per point per chron to determine accommodation. Tests with mcIterationtest indicate 20 is a reasonable minimum
    params.xcoVect = 0:150; % number of grid cells on the horizontal model profile
    params.gridCellLength = 1000; % Length of each grid cell in meters

    % Define the initial topography, including concave-up proximal fluvial profile, assuming a 151 point vector
    params.initialElevationVect = 50:-1:-100;
    modVect = ones(1,51); % Vector to adjust initial water depth profile, for example to produce exponential upper profile geometry
    modVect(1) = 10; % Most proximal wd multiplier
    for j =2:51 % Assume WD profile set up with shoreline at x index=51
        nextPoint = modVect(j-1) * 0.95; % Exponential decay rate
        if nextPoint > 1 % only include magnification values > 1 so distal profile maintains original linear gradient geometry
            modVect(j) = nextPoint;
        end
    end
    params.initialElevationVect(1:51) = params.initialElevationVect(1:51) .* modVect; % Apply the modficiation vector to the fluvial profile section
    params.initialElevationVect(131:151) = params.initialElevationVect(130); % Create a flat section at the distal end of the profile

    % % % Add shelf-break slope from x=60
    % % for j =61:80
    % %     params.initialElevationVect(j) = params.initialElevationVect(j-1)-20;
    % % end
    % % for j =81:151
    % %     params.initialElevationVect(j) = params.initialElevationVect(j-1)-1;
    % % end

     % Subsidence rates across the profile in m per My
    params.subsidRateVect = -100:2:200; % Simple rotational subsidence across whole model profile, zero subsidence node at x=50
    params.subsidRateVect(131:151) = params.subsidRateVect(130); % Create a constant subsidence section to keep the distal section flat

    % Sediment supply/production parameters
    params.maxCarbProductionRate = zeros(1,params.totalChrons);
    params.maxCarbProductionRate(1) = 5000; % Maximum carbonate production rate, as a vertical rate m per My
    productionRateIncrementAmplitude = 50; % 0 for simplest model, or 50 for interesting sed supply forcing;
    for t = 2:params.totalChrons
        params.maxCarbProductionRate(t)  = params.maxCarbProductionRate(t-1) + ((rand * productionRateIncrementAmplitude) - (productionRateIncrementAmplitude/2));
    end
    params.maxCarbProductionRate = max(params.maxCarbProductionRate, 0); % Replace any negative values with zeros
    fprintf("Maximum carbonate production rate from %2.1f to %2.1f m per My\n", min(params.maxCarbProductionRate), max(params.maxCarbProductionRate))

    % Depositional model parameters
    params.surfaceLight = 300; % from Bosscher and Schlager 1992
    params.extinctionCoeff = 0.05;
    params.saturatingLight = 2000;

    % Eustatic parameters e.g. long-term random walk and shorter-term periodic element
    params.eustaticCurve = zeros(1,params.totalChrons); % Initialise the eustatic curve paramter as a vector
    % Calculate the sinusoidal curve element
    sinCurveTime = 1:params.totalChrons * 2; % Sine curve needs to be longer than total model duration because it will be resampled to make asymmetric
    eustasyPeriod = 90; % zero-to-zero duration, in chrons
    eustasyAmplitude = 20; % zero to peak oscillation size, in m 20m for SL forced, zero for unforced
    sinCurve = (sin((sinCurveTime/eustasyPeriod) * 2 * pi) .* eustasyAmplitude);
    % Make the sea level curve asymmetric, with faster rise than fall
    asymmSeaLevel = sinCurve; % Start with the symmetrical sin curve
    counter = 0; % Counter variable to control the resampling of the sine curve to make it asymmetric
    asymmResampleCount = 5; % eg 5, removes 5 out of every 6 values of rising sea-level, so reduces period of rising limb by factor of 5/6
    for t = 1:numel(sinCurveTime)-1
        if sinCurve(t+1) - sinCurve(t) > 0 % Test if the t point on sea level curve is in the rising limb, so rate of change > 0
            if counter < asymmResampleCount % if counter indicates not yet removed target number of rising limb points...
                asymmSeaLevel(t) = NaN; % Set rising limb point as NaN so it will be removed from the curve
                counter = counter + 1;
            else
                counter = 0; % Reached the target number of points removed, current iteration point has been preserved, reset counter so next point will be removed
            end
        end
    end
    asymmSeaLevel = asymmSeaLevel(~isnan(asymmSeaLevel)); % Recalculate the sea level curve removing all the NaN values set in the rising limb sections
    % Calculate a random walk curve element
    randomWalkCurve = zeros(1,params.totalChrons);
    eustaticIncrementAmplitude = 5; % Zero for simplest model, 5 for Sl-forced model
    for t = 2:params.totalChrons
        randomWalkCurve(t) = randomWalkCurve(t-1) + ((rand * eustaticIncrementAmplitude) - (eustaticIncrementAmplitude/2));
    end
    % Calculate the final curve as the combination of periodic and random elements
    params.eustaticCurve = asymmSeaLevel(1:params.totalChrons) + randomWalkCurve;

    params.maxDeposRate = 1000.0; % in m per My, based on time-dependent rates from Sadler (1981) - check unit conversion in code!
    params.marineDeposP = 0.8;
    params.terrestrialDeposP = 0.0; % 0 because carbonate model
    params.maxErosionPGradientTerrestrial = 50.0; % gradient, meters per km, above which erosion probability is 1
    params.maxErosionPGradientMarine = 0.36; % gradient, meters per meter, above which erosion probability is 1; 0.36 is ~20 degrees
    params.maxErosionRateTerrestrial = 0.001; % m per My
    params.maxErosionRateMarine = 10; % m per My

    params.colourMap = [0, 1.0, 0.0, 0.0; ... % Water depth less than 0m, so subaerial exposure, red       
                        -20, 167/255, 199/255, 231/255; ... % Water depth 20m to 0m so high-energy carbonate
                        -1000, 0.7, 0.7, 0.7]; % Everything else below 20m down to 100m max, carbonate mud, light grey

    save(modelName,"params");

    if plotFlag
        figure
        tiledlayout(4,2); % Create a plot with 4 rows (first value) and 2 columns (second value)
        nexttile(1,[4,1])
        EMTVector = (1:params.totalChrons) * params.chronInterval;
        plot(params.eustaticCurve, EMTVector);
        grid on
        xlabel("Eustatic sea level (m)")
        ylabel("Elapsed model time (My)")
        nexttile(2,[4,1])
        plot(params.maxCarbProductionRate, EMTVector, "Color",[0.7,0,0]);
        grid on
        xlabel("Carbonate production rate (mMy-1)")
        ylabel("Elapsed model time (My)")
    end
end
