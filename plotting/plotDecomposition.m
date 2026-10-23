function fig = plotDecomposition(idxOutcome, irfs, cellChannelEffects, channelNames, ylabelText)
    % `plotDecomposition` Plot a decomposition of an IRF into transmission channels.
    %
    %   `fig = plotDecomposition(idxOutcome, irfs, cellChannelEffects, channelNames)`
    %   plots the channel effects on the outcome variable `idxOutcome` as
    %   stacked bars and the total effect as a line.
    %
    %   ## Arguments
    %   - `idxOutcome` (integer): Index of the outcome variable.
    %   - `irfs` (3D array): IRFs to a single shock, of size
    %     `(nVariables x 1 x (maxHorizon + 1))`. The total effect is taken from
    %     these IRFs.
    %   - `cellChannelEffects` (cell array): Channel effects, one entry per
    %     channel, each of the same size as `irfs`, as returned by `transmission`.
    %   - `channelNames` (string array): Names of the channels, used in the legend.
    %   - `ylabelText` (string, optional): Label of the y-axis. Defaults to
    %     `'Effect Size'`.
    %
    %   ## Returns
    %   - `fig` (Figure): Handle of the created figure.
    %
    %   See also `plotCompareDecompositions`, `transmission`

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

