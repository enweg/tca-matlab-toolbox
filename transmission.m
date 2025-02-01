function effects = transmission(from, arr1, arr2, q, method)
    % `transmission` Compute the transmission effect given a transmission condition.
    %
    %   `effects = transmission(from, arr1, arr2, q, method)` calculates the 
    %   transmission effect for a given transmission condition `q` using either 
    %   the `BOmega` method or the `irf` method. If `BOmega` is used, then 
    %   transmission effects will be calculated using the systems form 
    %   $x = Bx + \Omega\varepsilon$
    %
    %   ## Arguments
    %   - `from` (integer): Index of the shock.
    %   - `arr1` (matrix): 
    %       - If `method = "BOmega"`, this must be `B` from the systems form.
    %       - If `method = "irf"`, this must be structural `irfs` (technically
    %         only those of the shock that is being investigates (`from`)). 
    %         Has to be a IRF matrix. See `toTransmissionIrfs` for more information.
    %   - `arr2` (matrix): 
    %       - If `method = "BOmega"`, this must be `Omega` from the systems form.
    %       - If `method = "irf"`, this must be `irfsOrtho` (Cholesky IRFs) following 
    %         the ordering of the transmission matrix.
    %   - `q` (Q): A transmission condition. See also `Q`.
    %   - `method` (string): Specifies the calculation method:
    %       - `"BOmega"` uses the systems form.
    %       - `"irf"` uses only IRFs and can thus be used with local projections.
    %
    %   ## Returns:
    %   - `effects` (vector): A vector where entry `i` corresponds to the transmission 
    %     effect on variable `x_i`. If `x_k` is the variable in the transmission 
    %     condition with the highest index, all entries in the returned vector with 
    %     index less than `k` are `NaN`, since interpretation of those results is nonsensical.
    %
    %   ## Example:
    %   ```
    %   k = 6;
    %   h = 3;
    %   s = "(x1 | x2) & !x3";
    %   cond = makeCondition(s);
    %
    %   B = randn(k*(h+1), k*(h+1));
    %   Omega = randn(k*(h+1), k*(h+1));
    %
    %   effect = transmission(1, B, Omega, cond, "BOmega");
    %
    %   irfs = randn(k, k, h+1);
    %   irfsOrtho = randn(k, k, h+1);
    %
    %   irfs = toTransmissionIrfs(irfs);
    %   irfsOrtho = toTransmissionIrfs(irfsOrtho);
    %
    %   effect = transmission(1, irfs, irfsOrtho, cond, "irf");
    %   ```
    %
    %   ## Notes:
    %   - If `method = "BOmega"`, the function applies `transmissionBOmega`.
    %   - If `method = "irf"`, the function applies `transmissionIrfs`.
    %
    %   See also `transmissionBOmega`, `transmissionIrfs`, `makeCondition`, `through`, `notThrough`
    [varAnd, varNot, multiplier] = getVarNumsAndMultiplier(q);
    if isequal(method, "BOmega")
        disp('INFO: using method BOmega to compute transmission effect.');
        % apply the BOmega method
        % arr1 = B, arr2 = Omega
        effects = transmissionBOmega(from, arr1, arr2, varAnd, varNot, multiplier);
    elseif isequal(method, "irf")
        disp('INFO: using method irfs to compute transmission effect.');
        % apply the irf method 
        % arr1 = irfs, arr2 = irfsOrtho
        effectsCell = cell(1, length(varAnd)); 
        parfor ii = 1:length(varAnd) % Use parfor for parallel execution
            vAnd = varAnd{ii};
            vNot = varNot{ii};
            m = multiplier(ii);
            effectsCell{ii} = transmissionIrfs(from, arr1, arr2, vAnd, vNot, m);
        end
        effects = effectsCell{1};
        for ii = 2:length(effectsCell)
            effects = effects + effectsCell{ii};
        end
    else
        error("Unrecognised method. Must be either BOmega or irf");
    end
end
