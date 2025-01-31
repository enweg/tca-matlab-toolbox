function [hasContradiction, contradictions] = checkContradiction(varAnd, varNot)
    % `checkContradiction` Check for contradictions of the form `x1 & ~x1` (`x1 & !x1`).
    %
    %   `[hasContradiction, contradictions] = checkContradiction(varAnd, varNot)` 
    %   checks whether there is a contradiction where a variable appears in both 
    %   `varAnd` and `varNot`, meaning it is simultaneously required and forbidden.
    %
    %   ## Arguments
    %   - `varAnd` (vector or cell array of vectors): AND variable numbers obtained 
    %     from `getVarNumsAndMultiplier`.
    %   - `varNot` (vector or cell array of vectors): NOT variable numbers obtained 
    %     from `getVarNumsAndMultiplier`.
    %
    %   ## Returns
    %   - `hasContradiction` (logical): True if any contradictions exist.
    %   - `contradictions` (vector of logicals): A vector indicating which elements 
    %     yield a contradiction.
    %
    %   ## Notes
    %   - This function is used in `removeContradictions` to eliminate contradicting terms,
    %     which helps speed up the simplification process by reducing the total 
    %     number of terms.
    %
    %   See also `removeContradictions`, `getVarNumsAndMultiplier`

    % for vectors of integers
    if isnumeric(varAnd) && isnumeric(varNot)
        contradictions = ismember(varAnd, varNot);
        hasContradiction = any(contradictions);
        return;
    end

    % for cell arrays of integer vectors
    if iscell(varAnd) && iscell(varNot)
        contradictions = zeros(size(varAnd));
        hasContradiction = false;
        for ii=1:numel(varAnd)
            [iiHasContradiction, iiContradictions] = checkContradiction(varAnd{ii}, varNot{ii});
            hasContradiction = hasContradiction | iiHasContradiction;
            contradictions(ii) = any(iiContradictions);
        end
    else
        error('Invalid input types. Both inputs should be either numeric vectors or cell arrays of numeric vectors.');
    end
end

