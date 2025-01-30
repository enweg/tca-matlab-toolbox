function qOut = collectTerms(q)
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

