function idx = findVariableIndex(data, variable)
    % `findVariableIndex` Find column index of a variable in a table.
    %
    %   `idx = findVariableIndex(data, variable)` returns the index
    %   of `variable` in the table `data`.
    %
    %   ## Arguments
    %   - `data` (table): Dataset containing variable columns.
    %   - `variable` (char or integer): Name or index of the variable.
    %
    %   ## Returns
    %   - `idx` (integer): Column index of the variable.
    if isscalar(variable)
        idx = variable;
        return;
    end

    if ~istable(data)
        error("`data` must be a table if `variable` is not a scalar")
    end

    varnames = data.Properties.VariableNames;
    idx = find(cellfun(@(c) isequal(c, variable), varnames), 1, 'first');
end
