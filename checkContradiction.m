function [hasContradiction, contradictions] = checkContradiction(varAnd, varNot)
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

