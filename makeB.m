function B = makeB(As, Sigma, order, maxHorizon)
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
