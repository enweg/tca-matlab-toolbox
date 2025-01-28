function xIdx = mapY2XInt(i, t, K, order)
    xIdx = K * t + find(order == i, 1);
end
