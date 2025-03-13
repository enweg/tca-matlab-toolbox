function fig = plotDecomposition(idxOutcome, irfs, cellChannelEffects, channelNames, ylabelText)
    total = vec(irfs(idxOutcome, 1, :));
    channelEffects = vec(cellChannelEffects{1}(idxOutcome, 1, :));
    for ii = 2:length(cellChannelEffects)
        channelEffects = [channelEffects vec(cellChannelEffects{ii}(idxOutcome, 1, :))];
    end

    % Validate inputs
    if size(channelEffects, 1) ~= length(total)
        error('Total and channelEffects must have the same number of rows.');
    end
    
    % Set default ylabel if not provided
    if nargin < 5
        ylabelText = 'Effect Size';
    end
    
    % Define horizon
    horizon = length(total) - 1;
    
    % Create figure
    fig = figure;
    hold on;
    
    % Plot stacked bar chart
    bar(0:horizon, channelEffects, 'stacked', 'EdgeColor', 'none');
    
    % Plot total line
    plot(0:horizon, total, "-*", "LineWidth", 3, 'Color', "black");
    
    hold off;
    
    % Labels and legend
    xlabel('Horizon', 'FontSize', 15);
    ylabel(ylabelText, 'FontSize', 15);
    
    % Construct legend with all names
    legendNames = [channelNames, "Total"];
    lgd = legend(legendNames, ...
        'FontSize', 15, ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off' ...
    );
    
    % Formatting axes
    ax = gca;
    ax.FontSize = 15;
    ax.Box = 'on';
    ax.LineWidth = 2;
    grid on;
end

