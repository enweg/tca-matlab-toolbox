function effects = transmissionIrfs(from, irfs, irfsOrtho, varAnd, varNot, multiplier)
    effects = calculateQAndOnly(from, irfs, irfsOrtho, varAnd, 1); 

    combs = combinations(varNot);
    for i = 1:length(combs)
        c = combs{i};
        v = unique([varAnd, c]);
        effects = effects + (-1)^length(c) * calculateQAndOnly(from, irfs, irfsOrtho, v, 1);
    end

    effects = multiplier * effects;
end
