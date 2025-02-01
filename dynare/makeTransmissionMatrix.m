function T = makeTransmissionMatrix(variableOrder, vars)
    order = defineOrder(variableOrder, vars);
    T = eye(size(variableOrder));
    T = T(order, :);
end
