function q = notThrough(idx, horizons, order)
    % `notThrough` Construct a transmission condition excluding specific variables from the channel.
    %
    %   `q = notThrough(idx, horizons, order)` creates a transmission condition `Q` 
    %   where paths cannot pass through specified variables at given time horizons. 
    %   The variable indices refer to their positions in the original dynamic system 
    %   (e.g., in the SVARMA model), before applying the transmission matrix.
    %
    %   ## Arguments
    %
    %   **For a single variable:**
    %
    %   - `idx` (integer): Index of the variable that paths cannot go through, 
    %     using its original index in the dynamic system (before transmission ordering).
    %   - `horizons` (vector of integers): Time horizons at which the paths cannot 
    %     pass through the variable.
    %   - `order` (vector of integers): Variable ordering determined by the 
    %     transmission matrix.
    %
    %   **For multiple variables:**
    %
    %   - `idx` (vector of integers): Indices of the variables that paths cannot 
    %     go through, using their original indices in the dynamic system.
    %   - `horizons` (cell array of vectors): Each element specifies the time 
    %     horizons for the corresponding variable in `idx`.
    %   - `order` (vector of integers): Variable ordering determined by the 
    %     transmission matrix.
    %
    %   ## Returns
    %   - `q` (Q): A transmission condition.
    %
    %   ## Notes
    %   - The resulting transmission condition can be used in `transmission` to 
    %     compute the transmission effect.
    %
    %   ## Example
    %   ```
    %   % Non-contemporaneous monetary policy channel (Section 5.1 in Wegner etal 2024)
    %   q = notThrough(1, [0], 1:4);
    %
    %   % Anticipation channel of government defense spending (Section 5.2 in Wegner etal 2024)
    %   q = notThrough(2, 0:20, 1:4);
    %   ```
    %
    %   See also `through`, `transmission`

    if length(idx) > 1 && iscell(horizons)
        % Case when idx is a vector
        qs = arrayfun(@(ii) notThrough(idx(ii), horizons{ii}, order), 1:length(idx), 'UniformOutput', false);
        q = qs{1};
        for i = 2:length(qs)
            q = q & qs{i};
        end
    else
        % Case when idx is a scalar
        s = strjoin(arrayfun(@(h) sprintf('!y_{%d, %d}', idx, h), horizons, 'UniformOutput', false), ' & ');
        q = makeConditionY(s, order);
    end
end
