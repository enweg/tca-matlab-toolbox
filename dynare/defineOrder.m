function order = defineOrder(vars, options_)
    varsOriginal = dynareCellArrayToVec(options_.varobs);
    order = zeros(length(vars))
    for ii = 1:length(vars)
        order(ii) = find(varsOriginal == vars(ii));
    end
end
