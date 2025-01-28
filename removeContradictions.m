function q = removeContradictions(q)

    % If user does not define otherwise, we remove contradictions
    if exist('REMOVECONTRADICTIONS') == 1 && ~REMOVECONTRADICTIONS 
        return;
    end

    [varAnd, varNot, ~] = getVarNumsAndMultiplier(q);
    [hasContradiction, contradictions] = checkContradiction(varAnd, varNot);

    if hasContradiction
        if all(contradictions) % all terms contain some contradiction
            q = Q('T', 0); % This will later result in a zero
            return;
        end
        q = Q(q.vars(~contradictions), q.multiplier(~contradictions));
    end
end
