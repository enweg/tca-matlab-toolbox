function idx = findVariableIndex(data, variable)
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
