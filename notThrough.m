function q = notThrough(idx, horizons, order)
    if length(idx) > 1 && iscell(horizons)
        % Case when idx is a vector
        qs = arrayfun(@(ii) notThrough(idx(ii), horizons{ii}, order), 1:length(idx), 'UniformOutput', false);
        q = qs{1};
        for i = 2:length(qs)
            q = q & qs{i};
        end
    else
        % Case when idx is a scalar
        s = strjoin(arrayfun(@(h) sprintf('!y_{%d, %d}', idx, h), horizons, 'UniformOutput', false), ' & ');
        q = makeConditionY(s, order);
    end
end
