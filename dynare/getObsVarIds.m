function obsVar=getObsVarIds(M_, options_, oo_, varnames)
    % `getObsVarIds` Get Dynare variable IDs for given variable names.
    %
    %   `obsVar = getObsVarIds(M_, options_, oo_, varnames)` returns a vector of 
    %   Dynare variable IDs corresponding to the given variable names.
    %
    %   ## Arguments
    %   - `M_` (Struct): Returned by Dynare, containing model definitions.
    %   - `options_` (Struct): Returned by Dynare, containing simulation options.
    %   - `oo_` (Struct): Returned by Dynare, containing simulation results.
    %   - `varnames` (Cell array): Cell array of strings, where each string 
    %      is a variable name as specified in the .mod file.
    %
    %   ## Return
    %   - `obsVar` (Vector): IDs corresponding to the given variable names. 
    %      These IDs are the internal Dynare indices for the specified 
    %      variables.

    ids = repelem(1, length(varnames));
    endoname = dynareCellArray2Vec(M_.endo_names);
    for i=1:length(varnames)
        name = varnames(i);
        ids(i) = find(endoname == name);
    end
     
    obsVar = oo_.dr.inv_order_var(ids);
end
