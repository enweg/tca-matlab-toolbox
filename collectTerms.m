function q_out = collectTerms(q)
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
    non_zero_mult = mult ~= 0;
    vars = vars(non_zero_mult);
    mult = mult(non_zero_mult);

    % Create a new Q object with filtered variables and multipliers
    q_out = Q(vars, mult);
end

