function q = through(idx, horizons, order)
    % `through` Construct a transmission condition enforcing paths through specific variables.
    %
    %   `q = through(idx, horizons, order)` creates a transmission condition `Q` 
    %   where paths must pass through the variables specified in `idx`.
    %   The variable indices refer to their positions in the original dynamic system 
    %   (e.g., in the SVARMA model), before applying the transmission matrix.
    %
    %   ## Arguments
    %
    %   **For a single variable:**
    %
    %   - `idx` (integer): Index of the variable that paths must go through, 
    %     using its original index in the dynamic system (before transmission ordering).
    %   - `horizons` (vector of integers): Time horizons at which the paths must 
    %     pass through the variable.
    %   - `order` (vector of integers): Variable ordering determined by the 
    %     transmission matrix.
    %
    %   **For multiple variables:**
    %
    %   - `idx` (vector of integers): Indices of the variables that paths must 
    %     go through, using their original indices in the dynamic system.
    %   - `horizons` (cell array of vectors or vector of integer): If a single 
    %     vector of integers is provided, then it will be applied to each `idx`. 
    %     Alternatively, a cell array of integer vectors can be provided in which 
    %     case each element in the cell array applies to the respective element
    %     in `idx`.
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
    %   % Contemporaneous channel (Section 5.1 in Wegner)
    %   contemporaneous_channel = through(1, [0], 1:4);
    %
    %   % Effect through the federal funds rate in the first two periods
    %   q = through(1, [0, 1], 1:4);
    %
    %   % Effect through both the federal funds rate and output gap
    %   q = through([1, 2], {[0, 1], [0, 1]}, 1:4);
    %
    %   % Adjusting for a re-ordered system where the output gap comes first
    %   q = through([1, 2], {[0, 1], [0, 1]}, [2, 1, 3, 4]);
    %   ```
    %
    %   See also `notThrough`, `transmission`


    if length(idx) > 1 && ~iscell(horizons)
        % Allowing users to provide multiple indices that should be shut off
        % for the same horizons without having to provide the horizons separately.
        horizonsCell = cell(1, length(idx));
        for i = 1:length(idx)
            horizonsCell{i} = horizons;
        end
        horizons = horizonsCell;
    end
    if length(idx) > 1 && iscell(horizons)
        % Case when idx is a vector
        qs = arrayfun(@(ii) through(idx(ii), horizons{ii}, order), 1:length(idx), 'UniformOutput', false);
        q = qs{1};
        for i = 2:length(qs)
            q = q & qs{i};
        end
    else
        % Case when idx is a scalar
        s = strjoin(arrayfun(@(h) sprintf('y_{%d, %d}', idx, h), horizons, 'UniformOutput', false), ' & ');
        q = makeConditionY(s, order);
    end
end
