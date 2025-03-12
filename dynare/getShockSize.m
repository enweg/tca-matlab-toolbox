function [shockSize, idx] = getShockSize(M_, shockName)
    % `getShockSize` Obtain the standard deviation and index of a specified shock.
    %
    %   `[shockSize, idx] = getShockSize(M_, shockName)` computes the standard
    %   deviation (size) of a specified shock in a DSGE model estimated using Dynare.
    %   It also returns the index of the shock within the internal shock vector.
    %
    %   ## Arguments
    %   - `M_` (struct): Returned by Dynare.
    %   - `shockName` (string): The name of the shock whose size and index are required.
    %
    %   ## Returns
    %   - `shockSize` (double): The standard deviation of the specified shock.
    %   - `idx` (integer): The index of the shock in the internal shock vector.
    %

    shocks = dynareCellArrayToVec(M_.exo_names);
    idx = find(shocks == shockName);
    shockSize = sqrt(M_.Sigma_e(idx, idx));
end
