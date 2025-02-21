function plotModelAnalysis(params, strata, modelName)

    plotAccommodationDepositionErosionRates(params, strata, modelName)
    plotAccommodationSupplyRatio(params, strata, modelName)
    % plotForesetSlopeAnalysis(params, strata, modelName)
    % plotTopographicAnalysis(params, strata, modelName)
end

function plotAccommodationDepositionErosionRates(params, strata, modelName)

    if params.siliciclasticModel
   
        totalAccommPerIteration = (sum(strata.totalAccommodation,2) * 1000) / params.chronInterval; % Convert to m2 per My, same as sediment supply input parameter
        totalThicknessVectPerChron = strata.chrons(2:params.totalChrons,:) - strata.chrons(1:params.totalChrons-1,:);
        totalDeposPerChronPerMy = (sum(totalThicknessVectPerChron, 2) * 1000) / params.chronInterval;
        fluvialThicknessVectPerChron = (strata.chrons(2:params.totalChrons,:) - strata.chrons(1:params.totalChrons-1,:)) .* (strata.wdClass(2:params.totalChrons,:) == 1);
        fluvialDeposPerChronPerMy = (sum(fluvialThicknessVectPerChron, 2) * 1000) / params.chronInterval;
        marineDeposPerChronPerMy = totalDeposPerChronPerMy - fluvialDeposPerChronPerMy;
        totalErodedPerChronPerMy = (strata.erosionRecord * 1000) / params.chronInterval;

        figure
        semilogy((2:params.totalChrons) * params.chronInterval, totalAccommPerIteration(2:params.totalChrons), "color", [0,0,0], "LineWidth",2.0)
        hold on
        line((2:params.totalChrons) * params.chronInterval, totalDeposPerChronPerMy(1:params.totalChrons-1), "color", [0,0,0], "LineWidth",2.0, "LineStyle",":")
        line((2:params.totalChrons) * params.chronInterval, fluvialDeposPerChronPerMy(1:params.totalChrons-1), "color", [0,0.7,0], "LineWidth",2.0)
        line((2:params.totalChrons) * params.chronInterval, marineDeposPerChronPerMy(1:params.totalChrons-1), "color", [0,0,0.9], "LineWidth",2.0)
        line((2:params.totalChrons) * params.chronInterval, totalErodedPerChronPerMy(1:params.totalChrons-1), "color", [0.7,0,0.5], "LineWidth",2.0)
     
        xlabel("Elapsed model time (My)")
        ylabel("Accommodation, supply, deposition (m2My-1)")
        grid on
        
        legend("Total accommodation","Total deposition", "Fluvial deposition", "Marine deposition", "Erosion")
        % title(modelName)
    end

    drawnow
end

function plotAccommodationSupplyRatio(params, strata, modelName)

    figure

    % Plot the AS profiles, every numberOfProfiles, and record the shoreline x position for each profile
    numberOfProfiles = 10;
    chronIncrement = round((params.totalChrons-1) / numberOfProfiles);
    shorelinePlot = zeros(round(params.totalChrons / chronIncrement), 2); % to store xco value and as ratio
    j = 1;
    for t=2:chronIncrement:params.totalChrons
        lineColour = [0,0,t/params.totalChrons]; % Gradational from black (oldest) to blue (youngest)
        semilogy(params.xcoVect, strata.AS(t,:),"Color",lineColour, "lineWidth", 2, "DisplayName",sprintf("t=%d",t)); % Record legend label as chron number
        hold on
        % Store jth profile  shoreline x postion (-1 to start at 0) and the AS value at that position
        shorelinePlot(j,:) = [strata.shorelineXPos(t)-1, strata.AS(t,strata.shorelineXPos(t))];
        if strata.AS(t,strata.shorelineXPos(t))<=0
            fprintf("For t=%4.3f shoreline at %d km A=%6.5f A/S=%8.7f\n", t, shorelinePlot(j,1), strata.totalAccommodation(t,j), shorelinePlot(j,2))
        end
        j = j + 1;
    end

    % Plot patches of colour to mark A>S and A<S zones on the plot
    xco = [min(xlim), min(xlim), max(xlim), max(xlim)];
    yco = [1,max(ylim), max(ylim), 1];
    patch(xco, yco, [0,0,1],"FaceAlpha",0.25,"edgecolor","none","HandleVisibility","off")
    yco = [1,min(ylim), min(ylim), 1];
    patch(xco, yco, [1,0,0],"FaceAlpha",0.25,"edgecolor","none","HandleVisibility","off")

    % Plot circle at the shoreline x postions on each AS profile
    scatter(shorelinePlot(:,1), shorelinePlot(:,2), 50, [0.9647, 0.3765, 0.1765], "filled","HandleVisibility","off")

    xlabel("Distance (km)")
    ylabel("A/S")
    grid on
    % title(modelName)
    legend
    drawnow
end

function plotTopographicAnalysis(params, strata, modelName)

    topogRecord = zeros(params.totalChrons, 9);
    for t = 1:params.totalChrons
        topogRecord(t,:) = [min(strata.topogGradientRecord(t,:)), mean(strata.topogGradientRecord(t,:)), max(strata.topogGradientRecord(t,:)), ...
                min(strata.topogGradientRecord(t, 1:strata.shorelineXPos(t))), mean(strata.topogGradientRecord(t, 1:strata.shorelineXPos(t))), max(strata.topogGradientRecord(t,1:strata.shorelineXPos(t))), ...
                min(strata.topogGradientRecord(t, strata.shorelineXPos(t):numel(params.xcoVect))), mean(strata.topogGradientRecord(t,strata.shorelineXPos(t):numel(params.xcoVect))), max(strata.topogGradientRecord(t,strata.shorelineXPos(t):numel(params.xcoVect)))];
    end

    figure
    plotRange = 2:params.totalChrons;
    plot(plotRange * params.chronInterval, topogRecord(plotRange, 1:3), "color", [1,0,0.8])
    hold on
    plot(plotRange * params.chronInterval,topogRecord(plotRange, 4:6), "color", [0,1,0])
    plot(plotRange * params.chronInterval,topogRecord(plotRange, 7:9), "color", [0,0,1])
    xlabel("Elapsed model time (My)")
    ylabel("Topog gradients (m/km)")
    grid on
    legend("All min","all mean","all max", "T min","T mean","T max", "M min","M mean","M max")
    titleStr = sprintf("Gradient analysis:%s", modelName);
    title(titleStr)
    drawnow
end

function plotForesetSlopeAnalysis(params, strata, modelName)

    chronsToPlot = 10;
    chronIncrement = params.totalChrons / chronsToPlot;

    figure
    titleStr = sprintf("Foreset depos rate analysis:%s", modelName);
    title(titleStr)
    slopeRecordToPlot = strata.slopeDeposThickRecord(1:chronIncrement:params.totalChrons, :)';
    slopeDeposRateMax = max(strata.slopeDeposThickRecord, [], 2);
    slopeDeposRateMean = mean(strata.slopeDeposThickRecord, 2);
    subplot(3,1,1)
    plot(slopeRecordToPlot, "LineWidth",2)
    xlabel("Distance from shoreline (km)")
    ylabel("Depos rate (m per timestep)")
    legend
    subplot(3,1,2)
    plot(slopeDeposRateMean(2:params.totalChrons), "LineWidth",2)
    xlabel("Time (iteration)")
    ylabel("Mean depos rate (m per timestep)")
    subplot(3,1,3)
    plot(slopeDeposRateMax(2:params.totalChrons), "LineWidth",2)
    xlabel("Time (iteration)")
    ylabel("Max depos rate (m per timestep)")
 

    figure
    titleStr = sprintf("Foreset gradient analysis:%s", modelName);
    title(titleStr)
    slopeRecordToPlot = strata.slopeDeposGradientRecord(1:chronIncrement:params.totalChrons, :)';
    plot(slopeRecordToPlot, "LineWidth",2)
    slopeGradientMax = max(strata.slopeDeposGradientRecord, [], 2:params.totalChrons);
    slopeGradientMean = mean(strata.slopeDeposGradientRecord, 2:params.totalChrons);
    subplot(3,1,1)
    plot(slopeRecordToPlot, "LineWidth",2)
    xlabel("Distance from shoreline (km)")
    ylabel("Gradient (m per timestep)")
    legend 
    subplot(3,1,2)
    plot(slopeGradientMean(2:params.totalChrons), "LineWidth",2)
    xlabel("Time (iteration)")
    ylabel("Mean gradient (m per m)")
    subplot(3,1,3)
    plot(slopeGradientMax(2:params.totalChrons), "LineWidth",2)
    xlabel("Time (iteration)")
    ylabel("Max gradient (m per m)")
   
    drawnow
end

% function plotModelAnalytics(params, strata, modelName)
% 
%     figure
%     subplot(2,1,1)
%     plot(params.xcoVect, strata.totalAccommodation(2,:), "blue")
%     hold on
%     plot(params.xcoVect, strata.chrons(1,:), "black")
%     plot(params.xcoVect, strata.chrons(params.totalChrons-1,:), "black")
%     ylabel("Elevation & accommodation (m)")
%     yyaxis right
%     plot(params.xcoVect, strata.erosionPProfile, "r.","linewidth",2)
%     grid on
%     xlabel("Distance (km)")
%     ylabel("Erosion probability")
%     title(modelName)
% 
%     subplot(2,1,2)
%     plot(params.xcoVect, rad2deg(atan(strata.topogGradientRecord)))
%     xlabel("Xco (km)")
%     ylabel("Surface gradient (degrees)")
%     grid on
% end