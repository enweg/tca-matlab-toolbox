function effects = transmission(from, arr1, arr2, q, method, order)
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
    %   - `order` (vector, optional): variable ordering as defined by the transmission matrix.
    %
    %   ## Returns:
    %   - `effects` (vector): A vector where entry `i` corresponds to the transmission
    %     effect on variable `x_i`. If `x_k` is the variable with the highest index
    %     that enters the (simplified) transmission condition through an AND, then
    %     all entries with index `i <= k` are `NaN`, since the outcome variable must
    %     be ordered after all variables that paths must go through. Variables that
    %     only enter through a NOT never lead to `NaN`. See also `maskUninterpretable`.
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
    %   - If `order` is provided, the returned effects will be a 3D array 
    %     of dimension (nVariables, 1, horizons) where the variables are in the 
    %     original ordering (before applying the transmission matrix).
    %   - If `order` is not provided, the returned effects will be a (nVariable*horizons, 1)
    %     dimensional vector following the ordering after applying the transmission 
    %     matrix. This is similar to the matrix obtained via $(I - B)^{-1}\Omega$
    %     using the matrices from the systems form.
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

    % In the theory we assume that the outcome variable is ordered after all
    % variables in the transmission condition. A user might not be aware of
    % this and may interpret the effect on a variable ordered before (or equal
    % to) a conditioning variable. For NOT this is not a problem: a path to the
    % outcome cannot go through variables ordered after it, so the effect is
    % still correct. For AND this is a problem: the path cannot have gone
    % through the AND variable before reaching the outcome. All effects up to
    % and including the highest AND variable are thus set to NaN. This is done
    % here, for the full condition, and not for individual terms, since the
    % irf method computes NOT literals via AND-only terms (inclusion-exclusion).
    effects = maskUninterpretable(effects, varAnd);

    if nargin == 6
        % The case when order has been defined
        T = permmatrix(order);  % applying transpose inverses order
        k = length(order);
        horizon = size(effects, 1) / k;
        orderInv = T' * (1:k)';
        effect_tensor = zeros(k, 1, horizon);
        for h = 1:horizon
            tmp = effects(((h-1)*k+1):(h*k));
            effect_tensor(:, :, h) = tmp(orderInv);
        end
        effects = effect_tensor;
    end
end
