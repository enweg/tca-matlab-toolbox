function effects = maskUninterpretable(effects, varAnd)
    % `maskUninterpretable` Set transmission effects that must not be interpreted to NaN.
    %
    %   `effects = maskUninterpretable(effects, varAnd)` sets the effects on all
    %   variables `x_i` with `i <= k` to `NaN`, where `x_k` is the variable with
    %   the highest index that enters any term of the transmission condition
    %   through an AND (a variable that paths must go through).
    %
    %   The theory assumes that the outcome variable is ordered after all
    %   variables that the transmission condition puts restrictions on. A user
    %   might not be aware of this requirement and may thus interpret the effect
    %   on a variable that is ordered before (or is equal to) one of the
    %   conditioning variables.
    %   - For NOT this is not a problem. A path to `x_i` cannot go through any
    %     variable ordered after `x_i`, so a NOT on such a variable does not
    %     restrict the paths to `x_i`, and the effect on `x_i` itself is
    %     correctly zero if `x_i` is excluded. NOT variables therefore never
    %     trigger masking.
    %   - For AND this is a problem. A path to `x_i` with `i <= k` cannot go
    %     through `x_k` before reaching the outcome, so the question is
    %     nonsensical. The computed numbers have no meaningful interpretation
    %     and can even differ between the `BOmega` and the `irf` method. They
    %     are therefore set to `NaN`.
    %
    %   The masking is applied to the effect of the full condition and not to
    %   the effect of individual terms, since the `irf` method rewrites NOT
    %   literals as AND literals using inclusion-exclusion; masking those
    %   internal terms would wrongly mask effects of NOT conditions.
    %
    %   ## Arguments
    %   - `effects` (vector): Transmission effects in the ordering of the
    %     transmission matrix.
    %   - `varAnd` (cell array of vectors): Each cell contains a vector of variable
    %     indices that are included via AND in the term. Can be obtained using
    %     `getVarNumsAndMultiplier`.
    %
    %   ## Returns
    %   - `effects` (vector): The effects with the uninterpretable entries set to `NaN`.
    %
    %   ## WARNING
    %   Internal function. Should not be called by users directly.
    %
    %   See also `transmission`, `getVarNumsAndMultiplier`
    vars = cat(2, varAnd{:});
    if isempty(vars)
        return;
    end
    effects(1:max(vars)) = NaN;
end
