function qOut = collectTerms(q)
    % `collectTerms` Collect and sum multipliers for identical Boolean terms.
    %
    %   `qOut = collectTerms(q)` collects all terms $Q(b)$ where the Boolean statement 
    %   $b$ is the same and sums their multipliers. The result is a transmission 
    %   condition where each term appears only once, but with possibly different 
    %   multipliers (not restricted to ±1).
    %
    %   ## Arguments
    %   - `q` (Q): A transmission condition. See also `Q`.
    %
    %   ## Returns
    %   - `qOut` (Q): A new transmission condition where identical terms have been 
    %     combined with summed multipliers.
    %
    %   ## Example
    %   ```
    %   q = Q({"x1", "x1"}, [1, 1]);  
    %   collectTerms(q)
    %   % Output: Q({"x1"}, 2)
    %
    %   q = Q({"x1", "T", "x1"}, [1, 1, -1]);  
    %   collectTerms(q)
    %   % Output: Q({"T"}, 1)
    %   ```


    % Initialize a map to store terms
    terms = containers.Map();

    % Loop through variables and multipliers
    for i = 1:length(q.vars)
        v = q.vars{i};
        m = q.multiplier(i);

        if isKey(terms, v)
            terms(v) = terms(v) + m;
        else
            terms(v) = m;
        end
    end

    % Extract keys (variables) and values (multipliers)
    vars = keys(terms);
    mult = cell2mat(values(terms));

    % Filter out zero multipliers
    nonZeroMult = mult ~= 0;
    vars = vars(nonZeroMult);
    mult = mult(nonZeroMult);

    % if all mults are zero return zero
    if all(mult == 0)
        qOut = Q('T', 0);
        return;
    end

    % Create a new Q object with filtered variables and multipliers
    qOut = Q(vars, mult);
end

