function Omega = makeOmega(Phi0, Psis, Sigma, order, maxHorizon)
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
