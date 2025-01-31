function effects = transmissionIrfs(from, irfs, irfsOrtho, varAnd, varNot, multiplier)
    % transmissionIrfs Compute the transmission effect using the IRF method.
    %
    %   effects = transmissionIrfs(from, irfs, irfsOrtho, varAnd, varNot, multiplier) 
    %   calculates the transmission effect using impulse response functions (IRFs).
    %
    %   Arguments:
    %   - from (integer): Index of the shock.
    %   - irfs (matrix): Structural impulse responses (only the one of the `from`
    %     shock is needed). Must be a matrix which can be obtained from toTransmissionIrfs.
    %   - irfsOrtho (matrix): Cholesky IRFs that must use the ordering defined
    %     in the transmission matrix. 
    %   - varAnd (vector of integers): Indices of variables that must be included 
    %     (AND conditions). Can be obtained from getVarNumsAndMultiplier.
    %   - varNot (vector of integers): Indices of variables that must be excluded 
    %     (NOT conditions). Can be obtained from getVarNumsAndMultiplier.
    %   - multiplier (number): Multiplier associated with each term. Can be 
    %     obtained from getVarNumsAndMultiplier.
    %
    %   Returns:
    %   - effects (vector): A vector where entry `i` corresponds to the transmission 
    %     effect on variable `x_i`. 
    %
    %   Example:
    %   k = 6;
    %   h = 3;
    %   s = "x1";
    %   cond = makeCondition(s);
    %
    %   irfs = randn(k, k, h+1);
    %   irfsOrtho = randn(k, k, h+1);
    %
    %   irfs = toTransmissionIrfs(irfs);
    %   irfsOrtho = toTransmissionIrfs(irfsOrtho);
    %
    %   [varAnd, varNot, multiplier] = getVarNumsAndMultiplier(cond);
    %   effect = transmissionIrfs(1, irfs, irfsOrtho, varAnd{1}, varNot{1}, multiplier);
    %
    %   WARNING: 
    %  Internal function. Should not be called by users directly. 
    %
    %   See also transmission, getVarNumsAndMultiplier, makeCondition, through, notThrough
    effects = calculateQAndOnly(from, irfs, irfsOrtho, varAnd, 1); 

    combs = combinations(varNot);
    for i = 1:length(combs)
        c = combs{i};
        v = unique([varAnd, c]);
        effects = effects + (-1)^length(c) * calculateQAndOnly(from, irfs, irfsOrtho, v, 1);
    end

    effects = multiplier * effects;
end
