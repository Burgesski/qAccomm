function plotStrata(params, strata, totalChronsToPlot, modelName, animationCounter, lithoColourFlag, shorelineInc)

    fprintf("Plotting cross section and chronostrat diagram for %d chrons for model %s ...", totalChronsToPlot, modelName)

    if animationCounter == 0
        scrSz = get(0,'ScreenSize');
        figure('position',[1,1,scrSz(3)*0.9,scrSz(4)*0.9])
    else
        clf
    end

    subplot(2,1,1)
    set(gca,"fontsize", 12)

    yyaxis left
    depoThreshold = 0.01; % Only plot strata polygons if thickness > 1cm
    thickness = zeros(totalChronsToPlot, numel(params.xcoVect)); % Make sure thickness matrix is declared same size as strata chrons to simplify indexing
    thickness(2:totalChronsToPlot, :) = strata.chrons(2:totalChronsToPlot, :) - strata.chrons(1:totalChronsToPlot-1, :);
    thickness > depoThreshold;
    % Now draw polygons showing strata according to water depth of deposition
    for t=2:totalChronsToPlot
        for x=1:numel(params.xcoVect)-1

            if ~thickness(t,x) && thickness(t,x+1) % Wedge zero at x, thickening to the right
                xco = [params.xcoVect(x), params.xcoVect(x+1), params.xcoVect(x+1)];
                zco = [strata.chrons(t-1,x), strata.chrons(t-1,x+1),  strata.chrons(t,x+1)];
                if lithoColourFlag
                    patch(xco, zco, [0.5+(strata.sandProportion(t,x+1)*0.5), 0.5+(strata.sandProportion(t,x+1)*0.5), 0.5], "EdgeColor","none");
                else
                    patch(xco, zco, strata.wdColour(t,x+1,:),'EdgeColor','none');
                end
            end

            if thickness(t,x) && ~thickness(t,x+1) % Wedge thinning to zero at x+1
                xco = [params.xcoVect(x), params.xcoVect(x), params.xcoVect(x+1)];
                zco = [strata.chrons(t-1,x), strata.chrons(t,x),  strata.chrons(t,x+1)];
                if lithoColourFlag
                    patch(xco, zco, [0.5+(strata.sandProportion(t,x)*0.5), 0.5+(strata.sandProportion(t,x)*0.5), 0.5],"EdgeColor","none");
                else
                    patch(xco, zco, strata.wdColour(t,x,:),'EdgeColor','none');
                end
            end

            if thickness(t,x) && thickness(t,x+1) % Same lithology at two points
                xco = [params.xcoVect(x), params.xcoVect(x), params.xcoVect(x+1), params.xcoVect(x+1)];
                zco = [strata.chrons(t-1,x), strata.chrons(t,x),  strata.chrons(t,x+1), strata.chrons(t-1,x+1)];
                if lithoColourFlag
                    patch(xco, zco, [0.5+(strata.sandProportion(t,x)*0.5), 0.5+(strata.sandProportion(t,x)*0.5), 0.5],"EdgeColor","none");
                else
                    patch(xco, zco, strata.wdColour(t,x,:),'EdgeColor','none');
                end
            end
        end
    end

    hold on
    chronInterval = round(params.totalChrons / 10);
    if chronInterval > 0 % Trap situation of too few chrons to plot time lines - chron interval will be zero in this case
        % timelines = strata.chrons(1:chronInterval:totalChronsToPlot, :);
        for t = 1:chronInterval:totalChronsToPlot
            line(params.xcoVect, strata.chrons(t, :), "Color", [0,0,0]) % plot the chrons
        end
        line(params.xcoVect, strata.chrons(totalChronsToPlot, :), "Color", [0,0,0]) % plot the youngest chron also, to make sure final deposition is visible
    end

    % Plot sea-level
    line([strata.shorelineXPos(totalChronsToPlot), numel(params.xcoVect)], [params.eustaticCurve(totalChronsToPlot), params.eustaticCurve(totalChronsToPlot)], "color", [0,0,1]);

    % plot the shoreline positions if it is a siliclastic model, but not if it is a carbonate model
    if params.siliciclasticModel
        % Plot a line through every shoreline position
        shorelineElevationsLinIndex = sub2ind(size(strata.chrons), 1:totalChronsToPlot, strata.shorelineXPos(1:totalChronsToPlot));
        shorelineElevations = strata.chrons(shorelineElevationsLinIndex);
        shorelineXPosToPlot = strata.shorelineXPos(1:totalChronsToPlot) - 1; % Convert from index value to km value
        plot(shorelineXPosToPlot, shorelineElevations, "-", "color", [0.9647, 0.3765, 0.1765]);

        % plot a marker symbol every shorelineInc data points
        shorelineElevationsLinIndex = sub2ind(size(strata.chrons), 1:shorelineInc:totalChronsToPlot, strata.shorelineXPos(1:shorelineInc:totalChronsToPlot));
        shorelineElevations = strata.chrons(shorelineElevationsLinIndex);
        shorelineXPosToPlot = strata.shorelineXPos(1:shorelineInc:totalChronsToPlot) - 1; % Convert from index value to km value
        scatter(shorelineXPosToPlot, shorelineElevations, 10, [0.9647, 0.3765, 0.1765]); % blood orange open circls mark shoreline 
    end

    grid on
    
    xlim([0,150])
    minimumElevationInPlot = min(min(strata.chrons(totalChronsToPlot, :)));
    ylim([minimumElevationInPlot,100]) % zoomed view options for 100 chron simple sequence animation
    
    xlabel("Distance (km)")
    ylabel("Elevation (m)")

    % plot the accommodation profile using second y axis, scale on right
    yyaxis right
    ycoProfile = [strata.totalAccommodation(totalChronsToPlot,:),0, 0]; % Extract accomodation values along x axis for one timestep, and add zero value end and back-to-start points
    ycoProfile(1)=-0.01; % Models with no proximal erosion don't plot properly without this,not sure why - needs proper debugging
    xcoProfile = [params.xcoVect, max(params.xcoVect),0]; % Xco points, from axis vector, plus end point and back-to-start zero point to match ycoProfile values, ensure neat plot
    patch(xcoProfile, ycoProfile, [0 .2 .8], "EdgeColor","none", "facealpha", 0.2)
    ylim([min(min(strata.totalAccommodation)), max(max(strata.totalAccommodation))]);
    ylabel("Final timestep accommodation (m per ky)")

    % if params.siliciclasticModel % plot the mean sand proportion if it is a siliclastic model, but not if it is a carbonate model
    %     meanSandProportion = zeros(1,numel(params.xcoVect));
    %     for x = 2:numel(params.xcoVect)
    %         sumSandProportion = sum(strata.sandProportion(2:totalChronsToPlot, x));
    %         noneHiatusChronCount = sum(thickness(2:totalChronsToPlot,x) > 0.001);
    %         meanSandProportion(x) = sumSandProportion / noneHiatusChronCount;
    %     end
    % 
    %     % meanSandProportion = mean(strata.sandProportion(2:totalChronsToPlot, :),1); % Calculate the mean of values in dimension 1 so mean through time for each position of x
    %     line(params.xcoVect, meanSandProportion, "color", [0.9647, 0.3765, 0.1765], "LineWidth", 2);
    % end
    ax = gca();
    ax.YAxis(1).Color = [0 0 0]; % Set both y-axis colours to black
    ax.YAxis(2).Color = [0 0 0];
    yyaxis left % Leave the plot with the axis set to allow zoom on the cross-section plot

    subplot(2,1,2) % Chronostratigraphic diagram and time-series curves

    chronoTimeVect = (1:totalChronsToPlot) * params.chronInterval;

    if lithoColourFlag
        sandColourMap(:,:,1) = 0.5 + (strata.sandProportion(1:totalChronsToPlot,:) * 0.5);
        sandColourMap(:,:,2) = 0.5 + (strata.sandProportion(1:totalChronsToPlot,:) * 0.5);
        sandColourMap(:,:,3) = ones(totalChronsToPlot, size(sandColourMap,2)) .* 0.5;
        chronoColour = sandColourMap .* (thickness > 0.001);
    else
        chronoColour = strata.wdColour .* (thickness > 0.001);
    end

    for x = 1:numel(params.xcoVect)
         for t = 1:totalChronsToPlot
               if sum(chronoColour(t,x,:)) == 0
                    chronoColour(t,x,:) = [1,1,1];
               end
        end
    end

    imagesc(params.xcoVect, chronoTimeVect, chronoColour)
    set(gca,'YDir','normal') % invert the plot axis so that the chronostrat plots correctly, oldest at bottom, younging upwards
    hold on
    if params.siliciclasticModel % plot the shoreline positions if it is a siliclastic model, but not if it is a carbonate model
        % shorelineTimes = chronoTimeVect(1:shorelineInc: numel(chronoTimeVect));
        % shorelineXPosToPlot = strata.shorelineXPos(1:shorelineInc:totalChronsToPlot) - 1; % Convert from index value to km value
        % scatter(shorelineXPosToPlot, shorelineTimes, 10, [0.9647, 0.3765, 0.1765]); % blood orange open circls mark shoreline
    
        shorelineXPosToPlot = strata.shorelineXPos(1:totalChronsToPlot) - 1;
        line(shorelineXPosToPlot, chronoTimeVect, "color", [0.9647, 0.3765, 0.1765],"linewidth",2); % blood orange open circls mark shoreline
    end

    % Remove the x-axis tick labels so as not to interfer with sea and supply plots
    set(gca,'Xticklabel',[])
    
    % plot the sealevel curve and grid lines
    sedSeaPlotWidth = 10;
    seaLevelDataScaled = (max(params.xcoVect) - sedSeaPlotWidth) + ((params.eustaticCurve(1:totalChronsToPlot) - min(params.eustaticCurve)) / (max(params.eustaticCurve)-min(params.eustaticCurve)) * sedSeaPlotWidth);
    plot(seaLevelDataScaled, chronoTimeVect,"b");
    drawSeaSupplyCurveGridLines(params.eustaticCurve, max(params.xcoVect) - sedSeaPlotWidth, sedSeaPlotWidth, totalChronsToPlot * params.chronInterval, 5, 1, 1) % 3 ticks, 1 for rounded, 1 for text labels

    % Plot the sediment supply curve and grid lines
    totalDeposVolume = strata.totalDeposited * 100000;
    % Don't plot totalDeposScaled because plots almost same as total external supply for 1000 iteration plots
    % totalDeposScaled = ((totalDeposVolume(1:totalChronsToPlot) - min(totalDeposVolume)) / max(totalDeposVolume)) * sedSeaPlotWidth;
    % line(totalDeposScaled(2:totalChronsToPlot), (2:totalChronsToPlot) * params.chronInterval,"color",[0.0,0,0]);
    if params.carbonateModel
        maxCarbProdDataScaled  = (params.maxCarbProductionRate(1:totalChronsToPlot) - min(params.maxCarbProductionRate)) / (max(params.maxCarbProductionRate) - min(params.maxCarbProductionRate))  * sedSeaPlotWidth;
        line(maxCarbProdDataScaled(2:totalChronsToPlot), chronoTimeVect(2:totalChronsToPlot),"color",[0,0,0.5]);
    end
    if params.siliciclasticModel 
        if min(params.sedimentSupply) ~= max(params.sedimentSupply) % Time-variable supply because min and max are different so scale curve accordingly
            externalSupplyScale = (params.sedimentSupply(1:totalChronsToPlot) - min(params.sedimentSupply)) / (max(params.sedimentSupply) - min(params.sedimentSupply))  * sedSeaPlotWidth;
            line(externalSupplyScale(2:totalChronsToPlot), chronoTimeVect(2:totalChronsToPlot),"color",[0.5,0,0]);
        else % Otherwise supply is constant, so scale at sedSeaPlotWidth to plot at max right position on x axis
            externalSupplyScale = (params.sedimentSupply(1:totalChronsToPlot) / max(params.sedimentSupply))  * sedSeaPlotWidth;
            line(externalSupplyScale(2:totalChronsToPlot), chronoTimeVect(2:totalChronsToPlot),"color",[0.5,0,0], "linewidth",1.5);
        end
        
    end
    drawSeaSupplyCurveGridLines(totalDeposVolume, min(params.xcoVect), sedSeaPlotWidth, totalChronsToPlot * params.chronInterval, 5, 0, 0) % 3 ticks, 0 for unrounded , 0 for no ticks - numbers too big to be neat
    
    grid on
    xlim([0, max(params.xcoVect)])
    ylim([0, params.totalChrons * params.chronInterval])
    ylabel("Elapsed model time (My)")
    set(gca,"fontsize", 12)
    drawnow

    fprintf("Done\n");

    if totalChronsToPlot == params.totalChrons
        screenShotFName = strcat("screenShots/", modelName, ".png");
        fprintf("Exporting cross section and chronostrat diagram for %d chrons for model %s to file %s ...", totalChronsToPlot, modelName, screenShotFName)
        exportgraphics(gcf,screenShotFName)
        fprintf("Done\n\n")
    end
end

function drawSeaSupplyCurveGridLines(curveData, plotStartX, plotXWidth, plotEndTime, tickCount, roundedFlag, textFlag)

    upperLimit = ceil(max(curveData) / 10) * 10;
    lowerLimit = floor(min(curveData) / 10) * 10;
    plotRange =  upperLimit - lowerLimit;
    unroundedTickSize = plotRange/(tickCount-1);
    x = ceil(log10(unroundedTickSize) - 1);
    pow10x = x ^ 10;
    if pow10x ~= 0 && roundedFlag % Only round if not zero, and rounding flag set to 1. Avoid rounding large numbers - does not work, needs debugging
        roundedTickRange = ceil(unroundedTickSize / pow10x) * pow10x;
    else
        roundedTickRange = unroundedTickSize;
    end

    x = lowerLimit;
    for j=1:tickCount
        xPplotPosition = plotStartX + (((x - lowerLimit) / plotRange) * plotXWidth);
        line([xPplotPosition,xPplotPosition], [0, plotEndTime], "color", [0.1,0.1,0.1], "linestyle","-.")
        if textFlag
            text(xPplotPosition, 0, sprintf("%d",x),"HorizontalAlignment", "center", "VerticalAlignment", "top");
        end
        x = x + roundedTickRange;
    end
end