function q = makeConditionY(strY, order)
    strX = mapY2X(strY, order);
    q = makeCondition(strX);
end
