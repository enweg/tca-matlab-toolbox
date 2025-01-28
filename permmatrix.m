function P = permmatrix(order)
    P = eye(length(order));
    P = P(order, :);
end
