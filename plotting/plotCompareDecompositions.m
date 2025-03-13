function fig = plotCompareDecomposition(...
        idxOutcome, ...
        irfs, ...
        cellChannelEffects1, ...
        cellChannelEffects2, ...
        channelNames, ...
        decompositionNames, ...
        ylabelText ...
    )

    total = vec(irfs(idxOutcome, 1, :));

    if length(cellChannelEffects1) ~= length(cellChannelEffects2)
        error("Decompositions have difference number of channels.")
    end

    channelEffects1 = vec(cellChannelEffects1{1}(idxOutcome, 1, :));
    channelEffects2 = vec(cellChannelEffects2{1}(idxOutcome, 1, :));
    for ii = 2:length(cellChannelEffects1)
        channelEffects1 = [channelEffects1 vec(cellChannelEffects1{ii}(idxOutcome, 1, :))];
        channelEffects2 = [channelEffects2 vec(cellChannelEffects2{ii}(idxOutcome, 1, :))];
    end

    if size(channelEffects1, 1) ~= length(total) || size(channelEffects2, 1) ~= length(total)
        error('Total and channelEffects1 must have the same number of rows.');
    end
    
    % Set default ylabel if not provided
    if nargin < 7
        ylabelText = 'Effect Size';
    end
    
    % Define horizon
    horizon = length(total) - 1;

    % constants
    barWidth = 0.7;
    barGap = 0.05;
    barOffset = barWidth / 4 + barGap / 2;
    barTransparancy = 0.5;

    % Create figure
    fig = figure;
    hold on;
    
    barDecomposition1 = bar((0:horizon) - barOffset, channelEffects1, 'stacked', 'EdgeColor', 'none', 'BarWidth', barWidth / 2);
    barDecomposition2 = bar((0:horizon) + barOffset, channelEffects2, 'stacked', 'EdgeColor', 'none', 'BarWidth', barWidth / 2);
    for ii = 1:length(barDecomposition1)
        barDecomposition2(ii).FaceColor = barDecomposition1(ii).FaceColor; 
        barDecomposition2(ii).FaceAlpha = barTransparancy;
    end
    
    lineTotal = plot(0:horizon, total, "-*", "LineWidth", 3, 'Color', "black");
    hold off;
    
    % Labels and legend
    xlabel('Horizon', 'FontSize', 15);
    ylabel(ylabelText, 'FontSize', 15);
    
    % Construct legend with all names
    legendNames = [channelNames, "Total"];
    lgd = legend([barDecomposition1, lineTotal], ...
        legendNames, ...
        'FontSize', 15, ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off' ...
    );

    footnote = {...
        sprintf("High-opacity bars (left) correspond to the %s decomposition.", decompositionNames(1)), ...
        sprintf("Low-opacity bars (right) correspond to the %s decomposition.", decompositionNames(2)) ...
    };
    annotation('textbox', [0.1, 0.05, 0.8, 0.05], 'String', footnote, ...
                'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
    
    % Formatting axes
    ax = gca;
    ax.FontSize = 15;
    ax.Box = 'on';
    ax.LineWidth = 2;
    grid on;
end

