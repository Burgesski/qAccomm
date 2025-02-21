function strata = qAccomm(modelName, params, animationCounter, plotOneModelFlag, lithoColourFlag)

    rng(42);
    fprintf("\n\n\nRunning model %s for %d chrons ... \n", modelName, params.totalChrons)
    strata = calculateStrata(params, modelName, animationCounter, plotOneModelFlag, lithoColourFlag);
end

function strata = calculateStrata(params, modelName, animationCounter, plotOneModelFlag, lithoColourFlag)

    if animationCounter > 0
        scrSz = get(0,'ScreenSize');
        figure('position',[1,1,scrSz(3)*0.5,scrSz(4)*0.5])
        animationSubCounter = 0;
    end

    if params.siliciclasticModel == 1 % Silicilastic model with two grain sizes, sand and mud
        strata.accommodationSand = zeros(params.totalChrons, numel(params.xcoVect));
        strata.accommodationMud = zeros(params.totalChrons, numel(params.xcoVect));
        strata.sandProportion = zeros(params.totalChrons, numel(params.xcoVect));
    end
    strata.totalAccommodation = zeros(params.totalChrons, numel(params.xcoVect)); % Applies to all models, siliciclastic or carbonate 
    strata.shorelineXPos = zeros(1,params.totalChrons);
    strata.chrons = zeros(params.totalChrons, numel(params.xcoVect));
    strata.wdColour = zeros(params.totalChrons, numel(params.xcoVect),3);
    strata.AS = zeros(params.totalChrons, numel(params.xcoVect));
    strata.totalDeposited = zeros(1, params.totalChrons);
    strata.topogGradientRecord = zeros(params.totalChrons, numel(params.xcoVect));
    strata.erosionRecord = zeros(params.totalChrons, 1);
    strata.erosionPProfile = zeros(1, numel(params.xcoVect));
    strata.chrons(1,:) = params.initialElevationVect;
    strata.shorelineXPos(1) = findShorelineXPos(strata.chrons(1,:), params.eustaticCurve(1));
    strata.slopeDeposThickRecord = zeros(params.totalChrons, 20);
    strata.slopeDeposGradientRecord = zeros(params.totalChrons, 20);

    for t = 2:params.totalChrons

        if params.carbonateModel
            fprintf("EMT %4.3f (%d) SL:%3.2f SS (thick):%1.0f ", t*params.chronInterval, t, params.eustaticCurve(t), (params.maxCarbProductionRate(t) * params.chronInterval))
        else
            fprintf("EMT %4.3f (%d) SL:%3.2f SS (thick):%1.0f ", t*params.chronInterval, t, params.eustaticCurve(t), ((params.sedimentSupply(t) * params.chronInterval) / params.gridCellLength))
        end
        
        % Apply subsidence rate, calculate slopes and relative sealevel at each profile point
        strata.chrons(1:t-1,:) = strata.chrons(1:t-1,:) - (params.subsidRateVect * params.chronInterval); % Subside all the previous chrons
        strata.topogGradientRecord(t,:) = calculateTopographicGradients(params, strata, t);
        SLRelativeElev = strata.chrons(t-1, :) - params.eustaticCurve(t); % Record chron elevation relative to curent sea-level at time of deposition for chron t
        
        % Calculate water depth class and colour across the profile
        for x=1:numel(params.xcoVect)
            for j=size(params.colourMap,1):-1:1
                if SLRelativeElev(x) >= params.colourMap(j,1) % Params.colourMap first element is the chron elevation above/at which colour should be applied
                    strata.wdColour(t,x,:) = params.colourMap(j,2:4); % Params.colourMap elements 2:4 are the colour RGB values
                    strata.wdClass(t,x) = j; % Assign a class integer value so e.g. flvuial = 1, shoreface = 2, etc. Used in e.g. chronostrat diagram plotting
                end
            end
        end
        shorelineXPos = findShorelineXPos(strata.chrons(t-1,:), params.eustaticCurve(t));
        % strata.topogRecord(t,:) = [min(topogGradients), mean(topogGradients), max(topogGradients), ...
            % min(topogGradients(1:shorelineXPos)), mean(topogGradients(1:shorelineXPos)), max(topogGradients(1:shorelineXPos)), ...
            % min(topogGradients(shorelineXPos:numel(params.xcoVect))), mean(topogGradients(shorelineXPos:numel(params.xcoVect))), max(topogGradients(shorelineXPos:numel(params.xcoVect)))];

        % Loop along the grid and calculate accommodation value at each x point
        % Accommodation may be positive, indicating potential for deposition, or negative, indicating erosion needs to be calculated
        for x=1:numel(params.xcoVect)

            if params.siliciclasticModel == 1 % Calculate accommodation for the two grain size siliciclastic model
                pValues = [params.terrestrialSandDeposP, params.marineSandDeposP];
                [strata.accommodationSand(t,x), strata.erosionPProfile(x), bothFluvialAndMarine] = calculateOneAccommodationProbabilityProfile(params, pValues, t, strata.topogGradientRecord(t,x), strata.chrons(t-1,x), ((x==54 && t==6) * plotOneModelFlag));
                
                pValues = [params.terrestrialMudDeposP, params.marineMudDeposP];
                [strata.accommodationMud(t,x), strata.erosionPProfile(x), ~] = calculateOneAccommodationProbabilityProfile(params, pValues, t, strata.topogGradientRecord(t,x), strata.chrons(t-1,x), ((x==59 && t==30) * plotOneModelFlag));
                
                strata.totalAccommodation(t,x) = strata.accommodationSand(t,x) + strata.accommodationMud(t,x);
            else % Calculate accommodation for the carbonate model
                pValues = [params.terrestrialDeposP, params.marineDeposP];
                [strata.totalAccommodation(t,x), strata.erosionPProfile(x), ~] = calculateOneAccommodationProbabilityProfile(params, pValues, t, strata.topogGradientRecord(t,x), strata.chrons(t-1,x), ((x==59 && t==30) * plotOneModelFlag));
            end

            if strata.totalAccommodation(t,x) < 0
                strata.erosionRecord(t) = strata.erosionRecord(t) + abs(strata.totalAccommodation(t,x));
            end
            
            if plotOneModelFlag && bothFluvialAndMarine
                fprintf("Both marine and fluvial deposition in section at x=%d chron=%d\n",x,t );
            end
        end
        fprintf("Accomm %1.0f Erod %1.0f ", sum(strata.totalAccommodation(t,:)), strata.erosionRecord(t));
        
        % Calculate the deposition and hence the new chron surface, either siliciclastic or carbonate, for this iteration
        if params.siliciclasticModel
            strata = depositionErosionOneChronSiliciclastic(params, strata, t);
        elseif params.carbonateModel
            strata = depositOneChronCarbonate(params, strata, t);
        else
            fprintf("No deposition because no type specified in input parameters\n")
        end
        
        % Find and record the shoreline position after all processes in current iteration have been calculated, to best reflect the depositional
        % profile and generate steady shoreline evolution record
        strata.shorelineXPos(t) = findShorelineXPos(strata.chrons(t,:), params.eustaticCurve(t));
        oneChronThickness = strata.chrons(t,:) - strata.chrons(t-1,:);
        strata.totalDeposited(t) = sum(oneChronThickness);
        fprintf("Depos %1.0f (topog Max %2.1f min %2.1f) ", strata.totalDeposited(t), max(strata.chrons(t,:)), min(strata.chrons(t,:)));

        fprintf("\n")
        if animationCounter

            if animationSubCounter == 0 % First frame of an animation
                animationFName = strcat("animations/", modelName, ".gif");
                fprintf("Writing animation to %s\n", animationFName)
                exportgraphics(gcf, animationFName); % Write animation first frame, overwriting any previous animation
            end

            animationSubCounter = animationSubCounter + 1;

            if (animationSubCounter > animationCounter) || (t == params.totalChrons)
                plotStrata(params, strata, t, modelName, animationCounter, lithoColourFlag);
                exportgraphics(gcf, animationFName, Append=true); % Write next animation frame, appended to any previous animation
                animationSubCounter = 1;
            end
        end
    end
end

function shorelineXPos = findShorelineXPos(topogProfile, seaLevel)
% Find the point on the x profile where topographic height is closest to sea-level. This defines the shoreline position

    x = numel(topogProfile); % Set start x value to the distal end of the model profile
    while x > 0 && topogProfile(x) < seaLevel % loop along profile from distal to proximal while topog below sea-level
            x = x - 1;
    end

    shorelineXPos = x; % Set shoreline to most proximal position that is below sea-level, as found in loop above
    if x>1 && abs(seaLevel - topogProfile(x-1)) < abs(seaLevel - topogProfile(x)) 
        shorelineXPos = x-1; % Adjust the shoreline position by one grid cell if the adjacent more proximal elevation is closer to sea-level, above or below
    end
end

function topogGradients = calculateTopographicGradients(params, strata, t)
% Calculate the gradients along the t-1 chron surface, to be input for the calculation of the t chron surface elevations/accumulation etc
% return value as dimensionless gradient m/m

    topogGradients = zeros(1,numel(params.xcoVect));

    if t>3
        tStart = t-3;
    else
        tStart = 1;
    end

    for x = 1:numel(params.xcoVect)
        if x < 3
            gradXcoVect = 1:5;
        elseif x > numel(params.xcoVect) - 3
                gradXcoVect = numel(params.xcoVect)-5:numel(params.xcoVect);
            else
                gradXcoVect = x-2:x+2;
        end

        % topogReliefsToAverage = abs(strata.chrons(tStart:t-1, gradXcoVect(1:4)) - strata.chrons(tStart:t-1, gradXcoVect(2:5)));
        topogReliefsToAverage = strata.chrons(tStart:t-1, gradXcoVect(1:4)) - strata.chrons(tStart:t-1, gradXcoVect(2:5));
        topogGradients(x) = sum(topogReliefsToAverage,"all") / numel(topogReliefsToAverage); 
    end

    topogGradients = topogGradients ./ params.gridCellLength; % Convert to dimensionless gradient value
    topogGradients(topogGradients < 0) = 0; % Set negative gradients to zero
end

function [accommodation, erosionP, bothFluvialAndMarine] = calculateOneAccommodationProbabilityProfile(params, pValues, time, topogGradient, initialDeposSurfaceElev, plotFlag)

    totalDepositedThickness = zeros(1, params.mcIterations); % Useful to record the thickness of strata in each mc iteration model
    profileHeights = [0:0.1:50]'; % Create column vector or profile heights
    deposProfileHeightsProb = [profileHeights, zeros(numel(profileHeights),1)]; % Frequency count recorded how often deposition occurs at each height in depositional profile
    maxTotalThicknessSoFar = 0;
    bothFluvialAndMarine = 0;
    
    for mcCount = 1:params.mcIterations
        [totalDepositedThickness(mcCount), layerThickness,  erosionP] = calculateOneVerticalSection(params, pValues, time, topogGradient, initialDeposSurfaceElev);
        deposProfileHeightsProb(:,2) = deposProfileHeightsProb(:,2) + (deposProfileHeightsProb(:,1) < totalDepositedThickness(mcCount));

        if totalDepositedThickness(mcCount) > maxTotalThicknessSoFar
            maxTotalThicknessSoFar = totalDepositedThickness(mcCount);
            maxLayerThicknessSoFar = layerThickness;
        end
    end

    deposProfileHeightsProb(:,2) = deposProfileHeightsProb(:,2) / params.mcIterations; % Convert frequency count to probability

    % Code to check for combined fluvial and marine deposition in the same chron - used to find a good location for output for accumulation plots
    totalStratIterations = params.chronInterval / params.AATimeStep;
    deposElev = initialDeposSurfaceElev;
    for t = 2:totalStratIterations
        if layerThickness(t) > 0
            deposElev = deposElev + layerThickness(t);
            % Check if both fluvial and marine strata deposited AND that
            % marine thickness > fluvial thickness
            if initialDeposSurfaceElev < params.eustaticCurve(time) && deposElev > params.eustaticCurve(time) && params.eustaticCurve(time)-initialDeposSurfaceElev>deposElev-params.eustaticCurve(time)
                bothFluvialAndMarine = 1; 
            end
        end
    end
    
    if plotFlag
        plotOneModelVerticalSection(params, maxLayerThicknessSoFar, initialDeposSurfaceElev, deposProfileHeightsProb, totalDepositedThickness)
    end

    % Calculate final accommodation value as mean of the MC thickness distribution
    accommodation = mean(totalDepositedThickness);

    if plotFlag 
        fprintf("Plotted distribution accommodation=%5.4f\n", accommodation)
    end
end

function [totalDeposThickness, layerThickness,  erosionProbability] = calculateOneVerticalSection(params, pValues, time, topogGradient, deposSurfaceElev)
 
    iteration = 2;
    totalIterations = params.chronInterval / params.AATimeStep; % Calculating accommodation for one chron, so need iterations defined by aaTimeStep and chronDuration
    oneDeposThick = params.maxDeposRate * params.AATimeStep * rand;
    layerThickness = zeros(totalIterations,1);
    initialDeposElev = deposSurfaceElev; % Record this before any erosion occurs so that this initial value can be used to calculate total final erosion
    totalErosionThickness = 0;

    if deposSurfaceElev < params.eustaticCurve(time) % Surface below sea level
        [erosionProbability, deposProbability] = calculateMarineProbabilitiesFromGradient(params, pValues(2), topogGradient);
        oneErosionThick = params.maxErosionRateMarine * params.chronInterval * topogGradient;
        deposStartsBelowSeaLevel = 1;
    else % Surface above sea level
        [erosionProbability, deposProbability] = calculateTerrestrialProbabilitiesFromGradient(params, pValues(1), topogGradient);
        oneErosionThick = params.maxErosionRateTerrestrial * params.chronInterval * topogGradient;
        deposStartsBelowSeaLevel = 0;
    end

    while iteration < totalIterations

        % Check that p values and rates are still correct as depos setting changes during aggradation
        if deposSurfaceElev > params.eustaticCurve(time) && deposStartsBelowSeaLevel % deposition started below sea level but is now above sea-level so need to update p values etc
            [erosionProbability, deposProbability] = calculateTerrestrialProbabilitiesFromGradient(params, pValues(1), topogGradient);
            oneErosionThick = params.maxErosionRateTerrestrial * params.chronInterval * topogGradient;
        end

        if deposProbability > rand   % Deposition if calculated deposition probability is > random number 0-1
            deposSurfaceElev = deposSurfaceElev + oneDeposThick;
            layerThickness(iteration) = oneDeposThick;
        end

        if erosionProbability > rand % Erosion if calculated deposition probability is > random number 0-1

            deposSurfaceElev = deposSurfaceElev - oneErosionThick;

            erosionLoop = iteration; % Start erosion loop from current iteration because some deposition might have happened this iteration
            erosionRecord = oneErosionThick;
            while erosionRecord > 0 && erosionLoop > 1 % loop through older layers eroding down to oneErosionThick deep
                if layerThickness(erosionLoop) > 0.0
                    if erosionRecord > layerThickness(erosionLoop)
                        erosionRecord = erosionRecord - layerThickness(erosionLoop);
                        totalErosionThickness = totalErosionThickness + layerThickness(erosionLoop);
                        layerThickness(erosionLoop) = 0.0;
                    else
                        layerThickness(erosionLoop) = layerThickness(erosionLoop) - erosionRecord;
                        totalErosionThickness = totalErosionThickness + erosionRecord;
                        erosionRecord = 0; % This means no more erosion in the model loop 
                    end
                end
                erosionLoop = erosionLoop - 1;
            end
        end

        iteration = iteration + 1;
    end

    if deposSurfaceElev < initialDeposElev
        totalDeposThickness = deposSurfaceElev - initialDeposElev; % Negative accommodation case, erosion below initial deposition surface
    else
        totalDeposThickness =  sum(layerThickness); % Positive accommodation case, deposition above initial deposition surface
    end
end

function [terrestrialErosP, terrestrialDeposP] = calculateTerrestrialProbabilitiesFromGradient(params, terrestrialDeposP, gradient)
% Calculate terrestrial erosion and deposition probabiities
% Method assumes a minimum probability params.terrestrialDeposP at gradient = 0 rising via sigmoidal function to p=1 at gradient values params.maxErosionPGradientTerrestrial 
% This is done with a sigmoidal curvve from x=-10 to x=10, so this x axis has to be mapped onto the required gradient range
    
    % Define a sigmoidal function, asymptotic to x=-10 and x=10, that will be scaled horizontally to calculate the probability-gradient relationship
    sigmoidXVector = -10:0.1:10;
    sigmoidXVectorRange = max(sigmoidXVector) - min(sigmoidXVector);

    % Map gradients to the appropriate x position on the scaled sigmoid curve. Horizontal curve axes are scaled such that P rises sigmoidally
    % from 0 at x=0 to 1 at x=params.maxErosionPGradientTerrestrial
    xcoSigmoidTerrestrial = min(sigmoidXVector) + ((gradient ./ params.maxErosionPGradientTerrestrial) * sigmoidXVectorRange); 
    
    % now use the appropriately scaled gradient xco values, now in the range -10 to +10 to calculate each required probability
    terrestrialDeposP = (1 / (1 + exp(xcoSigmoidTerrestrial))) * terrestrialDeposP;
    terrestrialErosP = 1 - (1 / (1 + exp(xcoSigmoidTerrestrial)));
end


function [marineErosP, marineDeposP] = calculateMarineProbabilitiesFromGradient(params, marineDeposP, gradient)
% Calculate marine erosion and deposition probabiities
% Method assumes a minimum probability params.marineDeposP at gradient = 0 rising via sigmoidal function to p=1 at gradient values params.maxErosionPGradientMarine
% This is done with a sigmoidal curvve from x=-10 to x=10, so this x axis has to be mapped onto the required gradient range
    
    % Define a sigmoidal function, asymptotic to x=-10 and x=10, that will be scaled horizontally to calculate the probability-gradient relationship
    sigmoidXVector = -10:0.1:10;
    sigmoidXVectorRange = max(sigmoidXVector) - min(sigmoidXVector);

    % Map gradients to the appropriate x position on the scaled sigmoid curve. Horizontal curve axes are scaled such that P rises sigmoidally
    % from 0 at x=0 to to 1 at x=params.maxErosionPGradientMarine
    xcoSigmoidMarine = min(sigmoidXVector) + ((gradient ./ params.maxErosionPGradientMarine) * sigmoidXVectorRange); 
    
    % now use the appropriately scaled gradient xco values, now in the range -10 to +10 to calculate each required probability
    marineDeposP = (1 / (1 + exp(xcoSigmoidMarine))) * marineDeposP;
    marineErosP = 1 - (1 / (1 + exp(xcoSigmoidMarine)));
end

function strata = depositionErosionOneChronSiliciclastic(params, strata, t)

    sandToDeposit = (params.sedimentSupply(t) * params.supplySandProportion * params.chronInterval) / params.gridCellLength; % Convert sand supply in m2My-1 to meters per chron interval
    mudToDeposit = (params.sedimentSupply(t) * (1 - params.supplySandProportion) * params.chronInterval) / params.gridCellLength; % Convert mud supply in m2My-1 to meters per chron interval
    strata.chrons(t,:) = strata.chrons(t-1, :); % Copy the previous layer elevation as starting condition for current iteration layer
    x = 1;
    while (sandToDeposit > 0.001 || mudToDeposit > 0.001) && x < numel(params.xcoVect) % on the grid, and either sand or mud volume still available
       
        strata.AS(t,x) = strata.totalAccommodation(t,x) / (sandToDeposit + mudToDeposit); % calculate A/S at x for this time step t

        if strata.totalAccommodation(t,x) < 0 % Negative accommodation, so erosion

            sandToDeposit = sandToDeposit - (strata.totalAccommodation(t,x)); % Add eroded to sed supply, note use of minus because accommodation value is negative
            erodedSurfaceElevation = strata.chrons(t,x) + strata.totalAccommodation(t,x); % Calculate erosion surface elevation, note + because accommodation is negative
            for erosT = t:-1:1 % Lopp to erode current and older chrons if above erosion surface
                if strata.chrons(erosT,x) > erodedSurfaceElevation
                    strata.chrons(erosT,x) = erodedSurfaceElevation;
                end
            end

            % Record slope erosion to analyse development of graded slope etc
            xSlope = x - strata.shorelineXPos(t-1); %Xslope is distance of current profile point x from the shorline
            if xSlope > 0 && xSlope < 21 % So measure erosion on the 20 points seaward of the shoreline, assuming this will be the slope
                strata.slopeRecord(t,xSlope) = strata.totalAccommodation(t,x);
            end

        else % Positive accommodation, so deposition

            totalThicknessToDeposit = 0; % Set depositional thicknesses to zero until if statements below indicate otherwise
            oneSandThickness = 0;
            oneMudThickness = 0;

            if strata.accommodationSand(t,x) > 0.001                
                if sandToDeposit >= strata.accommodationSand(t,x) % If sediment supply exceeds accommodation at grid point x...
                    oneSandThickness = strata.accommodationSand(t,x); % Calculate thickness assuming fill all the accommodation         
                else
                    oneSandThickness = sandToDeposit * params.depositionalProportion; % Deposit the specified proportion of the remaining sediment budget at x
                end
                sandToDeposit = sandToDeposit - oneSandThickness; % Update sediment budget
                totalThicknessToDeposit = totalThicknessToDeposit + oneSandThickness;
            end

            if strata.accommodationMud(t,x) > 0.001
                if mudToDeposit >= strata.accommodationMud(t,x) % If sediment supply exceeds accommodation at grid point x...
                    oneMudThickness = strata.accommodationMud(t,x); % Fill all the accommodation
                else
                    oneMudThickness = (mudToDeposit * params.depositionalProportion); % Deposit half the remaining sediment budget at x
                end
                mudToDeposit = mudToDeposit - oneMudThickness; % Update mud sediment budget
                totalThicknessToDeposit = totalThicknessToDeposit + oneMudThickness;
            end

            if totalThicknessToDeposit > 0.001
                strata.chrons(t,x) = strata.chrons(t, x) + totalThicknessToDeposit;
                strata.sandProportion(t,x) = oneSandThickness / totalThicknessToDeposit;

                % Record slope deposition to analyse development of graded slope etc
                xSlope = x - strata.shorelineXPos(t-1); %Xslope is distance of current profile point x from the shorline
                if xSlope > 0 && xSlope < 21 % So measure deposition on the 20 points seaward of the shoreline, assuming this will be the slope
                    strata.slopeDeposThickRecord(t,xSlope) = totalThicknessToDeposit;
                    strata.slopeDeposGradientRecord(t, xSlope) = strata.topogGradientRecord(t,xSlope);
                end
            end
        end
        
        x = x + 1;
    end
end

function strata = depositOneChronCarbonate(params, strata, t)

    strata.chrons(t,:) = strata.chrons(t-1, :); % Copy the previous layer elevation as starting condition for current iteration layer
    waterDepth = params.eustaticCurve(t) - strata.chrons(t,:);
    for x=1: numel(params.xcoVect)-1
       
        if strata.totalAccommodation(t,x) < 0 % So negative accommodation, so erosion
            erodedSurfaceElevation = strata.chrons(t,x) + strata.totalAccommodation(t,x); % Calculate erosion surface elevation, note + because accommodation is negative
            for erosT = t:-1:1
                if strata.chrons(erosT,x) > erodedSurfaceElevation
                    strata.chrons(erosT,x) = erodedSurfaceElevation;
                end
            end

        else % Positive accommodation, so deposition
             % Bosscher and Schlager production variation with water depth
             sedToDeposit = params.maxCarbProductionRate(t) * params.chronInterval * tanh((params.surfaceLight * exp(-params.extinctionCoeff * waterDepth(x))) / params.saturatingLight);
                  
             if sedToDeposit >= strata.totalAccommodation(t,x) % If sediment supply exceeds accommodation at grid point x...
                strata.chrons(t,x) = strata.chrons(t, x) + strata.totalAccommodation(t,x); % Fill all the accommodation
            else
                strata.chrons(t,x) = strata.chrons(t, x) + sedToDeposit; % Deposit all the calculated carbonate thickness
            end
        end
    end
end
