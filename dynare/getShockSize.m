function shockSize = getShockSize(shockName, M_)
    % `getShockSize` Obtain the standard deviation of a specified shock.
    %
    %   `shockSize = getShockSize(shockName, M_)` computes the standard
    %   deviation (size) of a specified shock in a DSGE model estimated using Dynare.
    %
    %   ## Arguments
    %   - `shockName` (string): The name of the shock whose size and index are required.
    %   - `M_` (struct): Returned by Dynare.
    %
    %   ## Returns
    %   - `shockSize` (double): The standard deviation of the specified shock.
    %

    idx = getShockIdx(shockName, M_);
    shockSize = sqrt(M_.Sigma_e(idx, idx));
end
