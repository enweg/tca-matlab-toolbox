function Omega = makeOmega(Phi0, Psis, Sigma, order, maxHorizon)
    % `makeOmega` Construct the Omega matrix in the system representation $x = Bx + \Omega\varepsilon$.
    %
    %   ## Arguments
    %   - `Phi0` (matrix): Impact matrix for the shocks.
    %   - `Psis` (cell array of matrices): MA terms for the dynamic model.
    %   - `Sigma` (matrix): Covariance matrix of the reduced form errors.
    %   - `order` (vector): Ordering of variables.
    %   - `maxHorizon` (integer): Maximum IRF horizon.
    %
    %   ## Returns
    %   - `Omega` (matrix): Part of the systems form. 
    %
    %   See also `makeB`, `makeSystemsForm`

    % 1. Creating the transmission matrix
    T = permmatrix(order);
    K = size(Sigma, 1);
    
    % 2. Cholesky decomposition
    [L, D] = makeLD(T * Sigma * T');
    
    % 3. Compute Qt and Psis
    Qt = L * T * Phi0;
    Psis = cellfun(@(Psi) D * L * T * Psi * Phi0, Psis, 'UniformOutput', false);
    
    % 4. Creating Omega
    rowBlock = [cell2mat(fliplr(Psis)) D * Qt];
    Omega = zeros(K * (maxHorizon + 1), K * (maxHorizon + 1));
    Omega = slideIn(Omega, rowBlock);
end
