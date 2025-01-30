function B = makeB(As, Sigma, order, maxHorizon)
    % makeB Construct the matrix B in the system representation x = Bx + Omega * ε.
    %
    % Arguments:
    %   - As (cell array of matrices): Autoregressive coefficient matrices.
    %   - Sigma (matrix): Covariance matrix of the shocks.
    %   - order (vector): Ordering of variables given by the transmission matrix.
    %   - maxHorizon (integer): Maximum IRF horizon.
    %
    % Returns:
    %   - B (matrix): Part of the sytems representation.
    %
    %   See also makeOmega, makeSystemsForm.

    % 1. Creating the transmission matrix
    T = permmatrix(order);
    As = cellfun(@(A) T * A * T', As, 'UniformOutput', false);
    
    % 2. Cholesky decomposition
    [L, D] = makeLD(T * Sigma * T');
    K = size(Sigma, 1);  % assuming as many shocks as variables
    As = cellfun(@(A) D * L * A, As, 'UniformOutput', false);  % this gives DQ'A_i^* in the paper
    
    % 3. Creating B
    rowBlock = [cell2mat(fliplr(As)) eye(K) - D * L];
    B = zeros(K * (maxHorizon + 1), K * (maxHorizon + 1));
    B = slideIn(B, rowBlock);
end
