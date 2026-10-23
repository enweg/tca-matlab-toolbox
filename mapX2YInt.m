function [i, t] = mapX2YInt(xi, order)
    % `mapX2YInt` Map the index of a systems form variable `x` to a variable and horizon of the dynamic system.
    %
    %   `[i, t] = mapX2YInt(xi, order)` returns the variable index `i` and the
    %   horizon `t` such that `x_{xi}` corresponds to `y_{i, t}`. This is the
    %   reverse of `mapY2XInt`.
    %
    %   ## Arguments
    %   - `xi` (integer): Index of the variable in the systems form.
    %   - `order` (vector): Ordering of variables given by the transmission matrix.
    %
    %   ## Returns
    %   - `i` (integer): Index of the variable in the original ordering.
    %   - `t` (integer): Horizon of the variable. The contemporaneous period is
    %     horizon `0`.
    %
    %   See also `mapY2XInt`, `mapX2Y`

    K = length(order);
    t = floor((xi - 1) / K);
    i = xi - t * K;  % under the current transmission matrix
    i = order(i);    % under the original order
end
