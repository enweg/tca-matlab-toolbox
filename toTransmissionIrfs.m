function irfs = toTransmissionIrfs(irfs)
    % `toTransmissionIrfs` Transform a standard 3D IRF array into a 2D IRF matrix.
    %
    %   `irfs = toTransmissionIrfs(irfs)` converts an impulse response function (IRF) 
    %   array of dimensions `(n_variables × n_shocks × n_horizons)` into a 2D matrix. 
    %   The first horizon in the input corresponds to horizon 0.
    %
    % Arguments:
    %   - `irfs` (3D array): An IRF array of size `(n_variables × n_shocks × n_horizons)`, 
    %     where:
    %       - `n_variables` is the number of variables,
    %       - `n_shocks` is the number of shocks,
    %       - `n_horizons` is the number of forecast horizons.
    %
    %   ## Returns
    %   - `irfs` (2D matrix): A transformed IRF matrix of size 
    %     `(n_variables * n_horizons) × (n_variables * n_horizons)`. This is 
    %     equivalent to computing $(I - B)^{-1}Q$ using the systems form.
    %
    %   ## Example
    %   ```
    %   irfs3D = rand(4, 2, 10); % Example 3D IRF array with 4 variables, 2 shocks, and 10 horizons
    %   irfs2D = toTransmissionIrfs(irfs3D);
    %   ```
    %
    %   ## Notes
    %   - The first `n_shocks` column of the returned matrix are simply the standard 
    %     slices of the IRF 3D-array stacked vertically. The next `n_shocks` columns
    %     follow the same principle, but with the first `n_variables` columns being
    %     zero because the shocks are time t shocks which cannot affect time t-1
    %     variables. 
    [k, ns, maxHorizon] = size(irfs);
    maxHorizon = maxHorizon - 1;
    
    % Reshape the 3D array into a 2D matrix by stacking along the third dimension
    irfs = reshape(permute(irfs, [1, 3, 2]), k * (maxHorizon + 1), []);
    
    % Construct the output matrix
    irfsNew = [];
    for h = 0:maxHorizon
        padding = zeros(k * h, ns);
        truncated = irfs(1:end - k * h, :);
        irfsNew = [irfsNew, [padding; truncated]];
    end
    
    irfs = irfsNew;
end
