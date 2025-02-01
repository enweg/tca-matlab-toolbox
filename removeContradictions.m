function q = removeContradictions(q)
    % `removeContradictions` Remove contradicting terms from a transmission condition.
    %
    %   `q = removeContradictions(q)` removes terms that contain contradictions of 
    %   the form `x_i & !x_i`, which always evaluate to false and contribute zero 
    %   to the transmission effect. Behaviour of the function can be changed by 
    %   setting `REMOVECONTRADICTIONS=false` locally.
    %
    %   Arguments:
    %   - `q` (Q): A transmission condition. See also `Q` and `makeCondition`.
    %
    %   Returns:
    %   - If `REMOVECONTRADICTIONS` is set to `false`, the input `q` is returned unchanged.
    %   - If `REMOVECONTRADICTIONS` is `true` or not set:
    %       1. If all terms are contradicting, `Q("T", 0)` is returned, which 
    %          represents a transmission effect of zero.
    %       2. If some terms are non-contradicting, only the non-contradicting terms 
    %          are retained in the output.
    %
    %   ## Example
    %   ```
    %   setRemoveContradictions(true); % on by default
    %   
    %   q = Q('x1', 1);
    %   q = removeContradictions(q); % Returns q unchanged (no contradictions).
    %
    %   q = Q('x1 & !x1', 1);
    %   q = removeContradictions(q); % Returns Q('T', 0).
    %
    %   q = Q({'x1 & !x1', 'x1 & x2'}, [1, 1]);
    %   q = removeContradictions(q); % Returns Q('x2 & x1', 1).
    %
    %   setRemoveContradictions(false);
    %
    %   q = Q('x1 & !x1', 1);
    %   q = removeContradictions(q); % Returns the original q.
    %   ```
    %
    %   ## Notes
    %   - If `REMOVECONTRADICTIONS` is not explicitly set, contradictions are removed by default.

    % If user does not define otherwise, we remove contradictions
    global REMOVECONTRADICTIONS
    if isempty(REMOVECONTRADICTIONS)
        REMOVECONTRADICTIONS = true;
    end
    if ~REMOVECONTRADICTIONS 
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
