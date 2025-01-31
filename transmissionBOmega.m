function effects = transmissionBOmega(from, B, Omega, varAnd, varNot, multiplier)
    % `transmissionBOmega` Compute the transmission effect using the `BOmega` method.
    %
    %   `effects = transmissionBOmega(from, B, Omega, varAnd, varNot, multiplier)` 
    %   calculates the transmission effect based on the systems form $x = Bx + \Omega\varepsilon$. 
    %
    %   ## Arguments
    %   - `from` (integer): Index of the shock.
    %   - `B` (matrix): See the systems form.
    %   - `Omega` (matrix): See the systems form.
    %   - `varAnd` (cell array of vectors): Each cell contains a vector of variable 
    %     indices that must be included (AND conditions). Can be obtained using 
    %     `getVarNumsAndMultiplier`.
    %   - `varNot` (cell array of vectors): Each cell contains a vector of variable 
    %     indices that must be excluded (NOT conditions). Can be obtained using 
    %     `getVarNumsAndMultiplier`.
    %   - `multiplier` (vector of numbers): Multipliers associated with each term.
    %     Can be obtained using `getVarNumsAndMultiplier`.
    %
    %   ## Returns
    %   - `effects` (vector): A vector where entry `i` corresponds to the transmission 
    %     effect on variable `x_i`. 
    %
    %   ## Example
    %   ```
    %   k = 6;
    %   h = 3;
    %   s = "(x1 | x2) & !x3";
    %   cond = makeCondition(s);
    %
    %   B = randn(k*(h+1), k*(h+1));
    %   Omega = randn(k*(h+1), k*(h+1));
    %
    %   [varAnd, varNot, multiplier] = getVarNumsAndMultiplier(cond);
    %   effect = transmissionBOmega(1, B, Omega, varAnd, varNot, multiplier);
    %   ```
    %
    %  ##  WARNING 
    %  Internal function. Should not be called by users directly. 
    %
    %   See also `transmission`, `applyAndToB`, `applyNotToB`, `getVarNumsAndMultiplier`
    effects = cell(1, length(varAnd));
    for ii = 1:length(varAnd)
        vAnd = varAnd{ii};
        vNot = varNot{ii};
        m = multiplier(ii);

        BTilde = B;
        OmegaTilde = Omega;
        
        for v = vAnd
            [BTilde, OmegaTilde] = applyAndToB(BTilde, OmegaTilde, from, v);
        end
        
        for v = vNot
            [BTilde, OmegaTilde] = applyNotToB(BTilde, OmegaTilde, from, v);
        end
        
        effects{ii} = (eye(size(BTilde)) - BTilde) \ OmegaTilde(:, from);
        effects{ii} = m * effects{ii};
        
        % if isempty(vAnd) && isempty(vNot)
        %     return;
        % end
        
        if ~isempty(vAnd)
            effects{ii}(1:max(vAnd)) = 0;
        end
    end

    effects = sum(cat(2, effects{:}), 2);
    if ~isempty(varAnd)
        effects(1:max(cat(2, varAnd{:}))) = 0;
    end

end

