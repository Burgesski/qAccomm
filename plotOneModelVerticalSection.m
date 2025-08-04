function plotOneModelVerticalSection(params, maxLayerThicknessSoFar, initialDeposSurfaceElev, deposProfileHeightsProb, totalDepositedThickness)

    plotTimeThicknessAccumulation(params, maxLayerThicknessSoFar, initialDeposSurfaceElev, deposProfileHeightsProb) 
    plotThicknessDistribution(totalDepositedThickness)

end

function plotTimeThicknessAccumulation(params, layerThickness, initialSurfaceElevation, deposProfileHeightsProb)
    
    totalStratIterations = params.chronInterval / params.AATimeStep;
    accumulatedThickness = zeros(1,totalStratIterations);
    faciesColour = zeros(totalStratIterations,3);
    deposElev = zeros(totalStratIterations,1);
    deposProbability = zeros(totalStratIterations,1);
    
    % Define facies colours and extract same-height probabilities from deposProfileHeightsProb prior to any plotting
    accumulatedThickness(1) = 0;
    deposElev(1) = initialSurfaceElevation;
    [~, deposProfileIndex] = min( abs( deposProfileHeightsProb(:,1) - initialSurfaceElevation)); 
    deposProbability(1) = deposProfileHeightsProb(deposProfileIndex,2);
    for t = 2:totalStratIterations
        accumulatedThickness(t) = accumulatedThickness(t-1) + layerThickness(t);
        deposElev(t) = initialSurfaceElevation + accumulatedThickness(t-1); % t-1 because we want the elevation at the time of deposition, not after deposition
        for j=size(params.colourMap,1):-1:1
            if deposElev(t) > params.colourMap(j,1) % Params.colourMap first element is the chron elevation above which colour should be applied
                faciesColour(t,:) = params.colourMap(j,2:4); % Params.colourMap elements 2:4 are the colour RGB values
            end
        end
        [~, deposProfileIndex] = min( abs( deposProfileHeightsProb(:,1) - (deposElev(t) - initialSurfaceElevation))); 
        deposProbability(t) = deposProfileHeightsProb(deposProfileIndex,2);
    end

    stratigraphicCompleteness = sum(layerThickness > 0) / totalStratIterations;
    fluvialStratigraphicCompleteness = sum(layerThickness > 0 & deposElev > 0) / totalStratIterations;
    fprintf("Plotted vertical section has stratigraphic completeness %5.4f (%5.4f for fluvial deposition)\n", stratigraphicCompleteness, fluvialStratigraphicCompleteness)

    figure("position",[50,50, 600, 600])
    tiledlayout(3,3); % Create a plot with 3 rows (first value) and 3 columns (second value)
    
    % Create a time-cumulative thickness plot with solid colour marking intervals of deposition
    nexttile(1,[3,2]); % Create a plot covering 3 rows and 2 columns in the 3x3 figure space
    for t = 2:totalStratIterations
        line([(t-1) * params.AATimeStep, t * params.AATimeStep], [initialSurfaceElevation + accumulatedThickness(t) - layerThickness(t), initialSurfaceElevation + accumulatedThickness(t)], "color", "[0,0,0]", "linewidth",2)
        if layerThickness(t) > 0.0
            xco = [(t-1) * params.AATimeStep, (t-1) * params.AATimeStep, t * params.AATimeStep, t * params.AATimeStep];
            yco = [initialSurfaceElevation, initialSurfaceElevation + accumulatedThickness(t) - layerThickness(t), initialSurfaceElevation + accumulatedThickness(t), initialSurfaceElevation];
            patch(xco, yco, faciesColour(t,:), 'LineStyle','none');
        end
    end
    grid on;
    xlabel('Geological time (My)');
    ylabel('Elevation/thickness (m)');

     % Create a vertical section plotting probabiity of deposition and water-depth facies colour coding
    nexttile(3,[3,1]); % Create a plot covering 3 rows and 1 column in the 3x3 figure space
    for t = 2: totalStratIterations
        if layerThickness(t) > 0
            xco = [0, 0, deposProbability(t), deposProbability(t-1)];
            layerElevationBase = initialSurfaceElevation + accumulatedThickness(t-1);
            layerElevationTop = initialSurfaceElevation + accumulatedThickness(t);
            yco = [layerElevationBase, layerElevationTop, layerElevationTop, layerElevationBase];
            patch(xco, yco, faciesColour(t,:), 'LineStyle','none');
        end
    end
    grid on;
    xlabel('Probability of accumulation');
    ylabel('Elevation/thickness (m)');
    drawnow

    screenShotFName = strcat("vertAccumulation.png");
    fprintf("Exporting vertical accumulation plots to file %s ...", screenShotFName)
    exportgraphics(gcf,screenShotFName)
    fprintf("Done\n\n")
end

function plotThicknessDistribution(totalDepositedThickness)

    figure
    histogram(totalDepositedThickness, 10, "Normalization","probability")
    xlabel("Total accumulated thickness (m)")
    ylabel("Relative frequency")
    grid on
    drawnow
    fprintf("Accommodation PDF mean is %5.4f\n", mean(totalDepositedThickness));
    screenShotFName = strcat("accommPDF.png");
    fprintf("Exporting accommodation PDF plot to file %s ...", screenShotFName)
    exportgraphics(gcf,screenShotFName)
    fprintf("Done\n\n")
end