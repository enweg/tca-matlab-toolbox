function [i, t] = mapX2YInt(xi, order)
    K = length(order);
    t = floor((xi - 1) / K);
    i = xi - t * K;  % under the current transmission matrix
    i = order(i);    % under the original order
end
