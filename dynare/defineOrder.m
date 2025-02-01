function order = defineOrder(variableOrder, vars)
    order = zeros(size(variableOrder))
    for i=1:length(variableOrder)
        order(i) = find(vars == variableOrder{i})
    end
end
