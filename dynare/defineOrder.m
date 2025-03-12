function order = defineOrder(vars, options_)
    % `defineOrder` Determine the ordering (the transmission matrix) of
    % observed variables in a DSGE model. 
    %
    %   `defineOrder(vars, options_)` returns an ordering vector `order` that
    %   maps the variables in `vars` to their corresponding positions in the
    %   list of observed variables of a DSGE model estimated using Dynare. It 
    %   therefore defines the transmission matrix and can be used in `makeB`, 
    %   `makeOmega`, `makeSystemsForm`, `makeConditionY`, `notThrough`, `through`, 
    %   and `transmission`. 
    %
    %   ## Arguments
    %   - `vars` (vector): A list of observed variable names.
    %   - `options_` (struct): Returned by Dynare.
    %
    %   ## Returns
    %   - `order` (vector): A vector containing the indices of `vars` in the
    %     original observed variable list.
    %
    %   See also `transmission`, `through`, `notThrough`, `makeSystemsForm`.
    varsOriginal = dynareCellArrayToVec(options_.varobs);
    order = zeros(length(vars))
    for ii = 1:length(vars)
        order(ii) = find(varsOriginal == vars(ii));
    end
end
